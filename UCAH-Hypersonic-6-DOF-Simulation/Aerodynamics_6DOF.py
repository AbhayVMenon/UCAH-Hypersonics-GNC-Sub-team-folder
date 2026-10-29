# Aero dynamics 
import numpy as np

# --- geometry, referenced to body cross-section ---
from Vehicle_6DOF import A_PLAN_RATIO, S_REF

def getAero(mach, alpha):

    #CL, CD = 0.0, 0.0
    CN = 2.0 * alpha + A_PLAN_RATIO * 0.7 * 1.2 * alpha * abs(alpha)

    m = max(mach, 1.2)                    # base drag term is supersonic-only
    CA = 0.06 + 0.02 + 0.5 / m**2

    CL = CN * np.cos(alpha) - CA * np.sin(alpha)
    CD = CN * np.sin(alpha) + CA * np.cos(alpha)
    return CL, CD