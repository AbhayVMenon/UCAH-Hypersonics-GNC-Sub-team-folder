## Hub for vehicle parameters. All parameters are in SI units.

# NOTE: X_CG, X_REF, C_BAR, inertias, and fin limits are placeholders
import numpy as np

MASS   = 120.0                    # kg
G0     = 9.80665                  # m/s^2
D_REF  = 0.2                      # body diameter, m
LENGTH = 2.0                      # m
S_REF  = np.pi / 4 * D_REF**2     # reference area, m^2
C_BAR  = D_REF                    # reference length for moments, m

A_PLAN_RATIO = 6.4    # planform area / S_REF, for a 2.0 x 0.4 m body
X_CG   = .8                     # CG, m from nose
X_REF  = 1.0                      # point Aero's C_pitch is measured about, m from nose

I_XX   = 0.6                      # roll inertia, kg*m^2
I_YY   = 40.0                     # pitch inertia, kg*m^2
I_ZZ   = 40.0                     # yaw inertia, kg*m^2

# fin specs
DELTA_MAX  = np.deg2rad(20)       # fin deflection limit, rad
DELTA_RATE = np.deg2rad(200)      # fin rate limit, rad/s

FIN_AREA   = .01                # m^2, area of one fin
FIN_CENTER  = 0.05
FIN_FROM_CG_N = .5
FIN_FROM_CG_E = (D_REF/2)/np.sqrt(2)   # hinge on the body surface at 45 deg roll: radius = D_REF/2
FIN_FROM_CG_D = (D_REF/2)/np.sqrt(2)
FIN_TAU     = 0.02                     # fin actuator time constant, s
FIN_THETA   = 45            # angle from body axis to fin plane, rad

##### DECK/KEEL changes

DECK_FIN_HEIGHT     = .2               # m^2, area of one fin
DECK_FIN_WIDTH      = .1               # m^2, area of one fin
DECK_FIN_AREA       = DECK_FIN_HEIGHT * DECK_FIN_WIDTH
DECK_FIN_CENTER     = DECK_FIN_HEIGHT/2
DECK_FIN_FROM_CG_N  = .5
DECK_FIN_FROM_CG_E  = .1
DECK_FIN_FROM_CG_D  = .1
DECK_THETA          = 45            # angle from body axis to fin plane, rad

KEEL_FIN_HEIGHT    = .2               # m^2, area of one fin
KEEL_FIN_WIDTH     = .1               # m^2, area of one fin
KEEL_FIN_AREA      = KEEL_FIN_HEIGHT * KEEL_FIN_WIDTH
KEEL_FIN_CENTER    = KEEL_FIN_HEIGHT/2
KEEL_FIN_FROM_CG_N = .5
KEEL_FIN_FROM_CG_E = .1
KEEL_FIN_FROM_CG_D = .1
KEEL_THETA         = 45            # angle from body axis to fin plane, rad
TRAJECTORY = 0

vehicleParams = {
    "mass": MASS,
    "g0": G0,
    "d_ref": D_REF,
    "length": LENGTH,
    "s_ref": S_REF,
    "c_bar": C_BAR,
    "a_plan_ratio": A_PLAN_RATIO,
    "x_cg": X_CG,
    "x_ref": X_REF,
    "i_xx": I_XX,
    "i_yy": I_YY,
    "i_zz": I_ZZ,
    "delta_max": DELTA_MAX,
    "delta_rate": DELTA_RATE,
    "fin_area": FIN_AREA,
    "fin_center": FIN_CENTER,
    "fin_from_cg_n": FIN_FROM_CG_N,
    "fin_from_cg_e": FIN_FROM_CG_E,
    "fin_from_cg_d": FIN_FROM_CG_D,
    "fin_theta": FIN_THETA,
    "fin_tau": FIN_TAU
}
vehicleParams["I"]    = np.diag([I_XX, I_YY, I_ZZ])        # inertia matrix about the CG, body axes
vehicleParams["Iinv"] = np.linalg.inv(vehicleParams["I"])  # inverted once here, not every RK4 stage
vehicleParams["x_cg"]  = 1   # m aft of the nose tip
vehicleParams["c_ref"] = .2  