import numpy as np
from Atmosphere_6DOF import getAtm
from Frame_Math import earth_to_body, body_to_earth, quat_rate
from Vehicle_6DOF import vehicleParams, MASS
from AeroTable import getAero
from ControlSUrfaceModeling import fin_loads

N_STATE = 14
POS  = slice(0, 3)    # x, y, z                         NED position [m]
VEL  = slice(3, 6)    # u, v, w                         body velocity [m/s]
QUAT = slice(6, 10)   # q0..q3                          body->NED, scalar first
RATE = slice(10, 13)  # roll_rate, pitch_rate, yaw_rate [rad/s]
I_DELTA = 13          # delta                           fin deflection [rad]
G_NED = np.array([0.0, 0.0, 9.80665])


def derivs(t, state, getAero, getControls):
    """
    Equations of motion.

    Returns (state_dot, outputs):
      state_dot -- np.array, same shape as state; the ONLY thing RK4 integrates
      outputs   -- dict of quantities computed along the way (Mach, qbar, nz, ...)
                   that are worth logging. Computed once here so the sim loop
                   never has to recompute (and possibly disagree with) the physics.
    """
    # --- Unpack by name ---
    x, y, z = state[POS]
    u, v, w = state[VEL]
    quat    = state[QUAT]
    rate    = state[RATE]
    delta   = state[I_DELTA]
    c_ref = vehicleParams["c_ref"]
    x_cg = vehicleParams["x_cg"]

    m    = vehicleParams["mass"]
    I    = vehicleParams["I"]        # 3x3 inertia matrix about the CG
    Iinv = vehicleParams["Iinv"]     # precomputed once in Vehicle_6DOF
    S    = vehicleParams["s_ref"]
    
    h = -z
    V = np.linalg.norm([u, v, w])
    T, P, rho, a = getAtm(h)            # local air
    qbar  = 0.5 * rho * V**2

    mach = V / a
    alpha = np.arctan2(w, u)                # body angle-of-attack
    beta  = np.arcsin(v / V)                # body sideslip angle
    
    cBody = getAero(mach, h, np.degrees(alpha))   # only if getAero's grid is in degrees
    
    
    # # --- FULL 6DOF DIAGNOSTIC SCHEDULE ---
    # # --- MAX RANGE GLIDE SCHEDULE ---
    # if t < 20.0:
    #     # Phase 1: Low-Drag Ascent
    #     # Keep fins flat to punch through the dense lower atmosphere without bleeding velocity.
    #     delta = np.deg2rad(np.array([0.0, 0.0, 0.0, 0.0]))    
    # elif t > 20.0 and t < 80:
        
    #     delta = np.deg2rad(np.array([60.0, 60.0, 0.0, 0.0]))
    # elif t > 80.0 and t < 140:
    #     delta = np.deg2rad(np.array([0.0, 0.0, 60.0, 60.0]))
    # else:
    #     delta = np.deg2rad(np.array([0.0, 0.0, 0.0, 0.0]))

    deltas = getControls(t, state)   # control surface deflections from guidance

    # fullVehicleAero = getAero(mach, h, deltas, np.degrees(alpha))   # full aero with fins


    F_N = qbar * S * cBody.CN               # perpendicular to the body axis, positive "up" (-z)
    F_A = qbar * S * cBody.CA               # along the body axis, positive aft (-x)

    # Calculate Fin loads with DAMPING (rate and V included)
    F_Fins, M_Fins = fin_loads(qbar, alpha, beta, delta, mach=mach, rate=rate, V=V)

    # Isolate the body aero forces to calculate the nose moment
    F_aero = np.array([-F_A, 0.0, -F_N])
    
    # --- Aero moment: table gives it about the nose, the EOMs need it about the CG ---
    M_nose = qbar * S * c_ref * np.array([0.0, cBody.Cm_nose, 0.0])
    r_nose_from_cg = np.array([x_cg, 0.0, 0.0])          # nose is x_cg ahead of the CG
    
    # Cross ONLY the body aero forces with the nose moment arm, then add fin moments
    M_body = M_nose + np.cross(r_nose_from_cg, F_aero) + M_Fins  

    # Combine all forces for translational acceleration (u_dot, v_dot, w_dot)
    F_body = F_aero + F_Fins   

    # --- 3-AXIS AERODYNAMIC DAMPING (Empirical Estimates) ---
    C_lp = -20.0   # Roll damping: mostly from fins resisting spin
    C_mq = -400.0  # Pitch damping: long body resisting rotation
    C_nr = -150.0  # Yaw damping: identical to pitch for symmetric bodies
    
    C_damp = np.array([C_lp, C_mq, C_nr])
    
    # Damping moment = qbar * S * d_ref * (d_ref / 2V) * C_damp * [p, q, r]
    damping_moment = qbar * S * c_ref * (c_ref / (2.0 * V)) * C_damp * rate
    
    M_body += damping_moment

    # EOMs
    vel = state[VEL]
    state_dot = np.zeros(N_STATE)
    state_dot[VEL]  = F_body / m - np.cross(rate, vel) + earth_to_body(G_NED, quat)
    state_dot[RATE] = Iinv @ (M_body - np.cross(rate, I @ rate))
    state_dot[POS]  = body_to_earth(vel, quat)
    state_dot[QUAT] = quat_rate(quat, rate)
    state_dot[I_DELTA] = 0.0

    outputs = {"h": h, "V": V, "mach": mach, "qbar": qbar,
                "alpha": alpha, "beta": beta,
                "nz": -F_body[2] / (m * G_NED[2])}   # normal load factor in g's (aero only; body z is down)

    return state_dot, outputs