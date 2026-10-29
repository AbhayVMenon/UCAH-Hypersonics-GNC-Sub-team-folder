# Sim loop
from AeroTable import getAero
from dynamics_6DOF import derivs, QUAT
from Guidance_6DOF import getControls
from RK4_6DOF import rk4
import numpy as np

STATE_KEYS  = ["x", "y", "z", "u", "v", "w", "q0", "q1", "q2", "q3",
               "roll_rate", "pitch_rate", "yaw_rate", "delta"]
OUTPUT_KEYS = ["h", "V", "mach", "qbar", "alpha", "beta", "nz"]      # pulled from derivs' outputs dict, in this order
COLUMNS     = ["t"] + STATE_KEYS + OUTPUT_KEYS              # history column order (display_6DOF unpacks in this order)


def run6DOF(initialState, vehicleParams=None, simParams=None):
   
    simParams = simParams or {}
    dt    = simParams.get("dt", 0.01)       # integration timestep, s
    t_max = simParams.get("t_max", 600.0)   # 10 min cutoff, s

    # bind aero + guidance so RK4 only sees f(t, state)
    f = lambda t, s: derivs(t, s, getAero, getControls)

    t = 0.0
    state = np.array(initialState, dtype=float)
    history = []                            # one row logged per timestep

    while state[2] < 0.0 and t < t_max:     # until ground impact or timeout
        if not np.all(np.isfinite(state)) or state[0] < 0.0:
            print(f"diverged at t={t:.2f}")
            break

        new_state, out = rk4(f, t, state, dt)   # out = derivs outputs at (t, state)
        new_state[QUAT] /= np.linalg.norm(new_state[QUAT])   # pull quaternion back to unit length

        # log the row BEFORE advancing, so state and outputs both match this t
        history.append([t, *state, *(out[k] for k in OUTPUT_KEYS)])

        state = new_state
        t += dt
    # this interpolates to find the true impact time and state, and logs one last row at that point 
    # previously, the state would propogate through the ground becaust of off timing
    if state[2] >= 0.0 and history:
        n = len(STATE_KEYS)
        prev_state = np.array(history[-1][1:1+n])   # last above-ground state
        h_prev, h_next = prev_state[2], state[2]

        frac = h_prev / (h_prev - h_next)          # 0..1: how far into the last step h hit 0
        t = (t - dt) + frac * dt                   # impact time
        state = prev_state + frac * (state - prev_state)   # impact state (h is now exactly 0)

        _, out = f(t, state)                       # one derivs call for consistent Mach, qbar, nz at impact
        history.append([t, *state, *(out[k] for k in OUTPUT_KEYS)])

    return np.array(history), state, t
