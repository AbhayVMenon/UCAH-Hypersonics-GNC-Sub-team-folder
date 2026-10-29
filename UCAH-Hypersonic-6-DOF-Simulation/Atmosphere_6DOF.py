import numpy as np
from scipy.interpolate import CubicSpline

# --- 1. PRE-COMPUTE THE BLENDED SPLINES ---
def _generate_exact_atmosphere(alt):
    """Original exact 1976 Standard Atmosphere logic to generate truth data."""
    H = 6356766.0 * alt / (6356766.0 + alt)

    if H < 11000.0:
        T = 288.15 - 0.0065 * H
        P = 101325.0 * (T / 288.15) ** 5.2559
    elif H < 20000.0:
        T = 216.65
        P = 22632.0 * np.exp(-9.80665 * (H - 11000.0) / (287.0528 * T))
    elif H < 32000.0:
        T = 216.65 + 0.0010 * (H - 20000.0)
        P = 5474.9 * (T / 216.65) ** -34.163
    elif H < 47000.0:
        T = 228.65 + 0.0028 * (H - 32000.0)
        P = 868.02 * (T / 228.65) ** -12.201
    elif H < 51000.0:
        T = 270.65
        P = 110.91 * np.exp(-9.80665 * (H - 47000.0) / (287.0528 * T))
    elif H < 71000.0:
        T = 270.65 - 0.0028 * (H - 51000.0)
        P = 66.939 * (T / 270.65) ** 12.201
    elif H < 84852.0:
        T = 214.65 - 0.0020 * (H - 71000.0)
        P = 3.9564 * (T / 214.65) ** 17.082
    else:
        T = 186.946
        P = 0.37338 * np.exp(-9.80665 * (H - 84852.0) / (287.0528 * T))
        
    return T, P

# Generate truth tables at 500m resolution from 0 to 100km
_alts = np.linspace(0, 100000, 201)  
_T_exact = np.zeros_like(_alts)
_lnP_exact = np.zeros_like(_alts)

for i, a in enumerate(_alts):
    t, p = _generate_exact_atmosphere(a)
    _T_exact[i] = t
    # Splining the natural log of pressure respects the exponential decay curve perfectly
    _lnP_exact[i] = np.log(p) 

# Create the C2 continuous interpolators
_T_spline = CubicSpline(_alts, _T_exact)
_lnP_spline = CubicSpline(_alts, _lnP_exact)

# --- 2. DROP-IN REPLACEMENT FUNCTION ---
def getAtm(alt):
    """
    Blended 1976 Standard Atmosphere (valid 0-100 km).
    Uses C2 continuous splines to prevent sharp derivative jumps in Mach/Speed of Sound.
    """
    # Evaluate smooth temperature and pressure from the splines
    T = _T_spline(alt)
    P = np.exp(_lnP_spline(alt))
    
    # Recompute dependent variables to maintain exact thermodynamic consistency
    rho = P / (287.0528 * T)          # ideal gas law
    a = np.sqrt(1.4 * 287.0528 * T)   # speed of sound
    
    # Ensure return types match standard python floats if a single float was passed,
    # while automatically supporting numpy array inputs if needed later.
    if np.isscalar(alt):
        return float(T), float(P), float(rho), float(a)
        
    return T, P, rho, a