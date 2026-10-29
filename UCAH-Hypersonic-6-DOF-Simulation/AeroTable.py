import os
from collections import namedtuple
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, CubicSpline

# Default table lives next to this file, so the lookup works from any working directory
DEFAULT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "HARVaerotable.csv")

# The three table inputs, in the order the grid is built
GRID_COLUMNS = ["Mach", "alt_m", "alpha_deg"]

# Columns we interpolate. q_kPa is left out: it's a flight condition, not an aero coefficient.
COEFF_COLUMNS = ["CN", "CA", "CY", "CL", "CD", "CDpres", "CDfric", "LD",
                 "Cm_nose", "Cl", "Cn", "xcp_m", "xcp_over_L",
                 "CYbeta", "Clbeta", "Cnbeta"]

# Named result, so you write c.CN instead of c[0]
AeroCoeffs = namedtuple("AeroCoeffs", COEFF_COLUMNS)

# Loaded tables and interpolators, keyed by csv path
_tables = {}
_interpolators = {}


def _loadTable(csv_path):
    """Read the CSV into breakpoint axes and a (Mach, alt, alpha, coeff) value grid."""
    df = pd.read_csv(csv_path)

    # The unique breakpoints along each axis
    axes = [np.sort(df[col].unique()).astype(float) for col in GRID_COLUMNS]
    shape = tuple(len(a) for a in axes)

    if len(df) != np.prod(shape) or df.duplicated(GRID_COLUMNS).any():
        raise ValueError("Aero table is not a complete Mach x alt x alpha grid")

    df = df.sort_values(GRID_COLUMNS)
    values = df[COEFF_COLUMNS].to_numpy(dtype=float).reshape(shape + (len(COEFF_COLUMNS),))
    
    # --- NEW ADDITION: Clean broken data before cubic interpolation ---
    if not np.isfinite(values).all():
        print(f"Warning: Missing (NaN) or infinite values found in {csv_path}. Safely replacing with 0.0.")
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    # ------------------------------------------------------------------

    # Pre-build SciPy RegularGridInterpolators for 3D Linear and Cubic methods
    interp_linear = RegularGridInterpolator(axes, values, method='linear', bounds_error=False, fill_value=None)
    interp_cubic = RegularGridInterpolator(axes, values, method='cubic', bounds_error=False, fill_value=None)
    
    return axes, values, interp_linear, interp_cubic


def getAero(mach, alt_m, alpha_deg, csv_path=DEFAULT_CSV, method='linear', warn_outside=True):
    """
    Return every aero coefficient at this flight condition.
    
    Available interpolation methods:
    - 'linear': 3D Linear interpolation (Standard for Cn, Ca, Cl, Cd, Cm)
    - 'cubic': 3D Cubic interpolation (Applies across M, H, alpha)
    - 'spline_alpha': Uses 2D cubic interpolation for (M, H) and a Cubic Smoothing Spline for alpha.
    """
    if csv_path not in _tables:
        _tables[csv_path] = _loadTable(csv_path)
        
    axes, values, interp_linear, interp_cubic = _tables[csv_path]

    point = np.array([mach, alt_m, alpha_deg])
    
    if warn_outside:
        for name, axis, x in zip(GRID_COLUMNS, axes, point):
            if x < axis[0] or x > axis[-1]:
                print(f"AeroTable: {name} = {x} outside table "
                      f"[{axis[0]}, {axis[-1]}], clamped to edge")
                
    # Clamp inputs to table edges to prevent wild extrapolation
    point = [np.clip(point[i], axes[i][0], axes[i][-1]) for i in range(3)]

    if method == 'linear':
        # 3D Linear Interpolation
        coeffs = interp_linear(point)[0]
        
    elif method == 'cubic':
        # 3D Cubic Interpolation
        coeffs = interp_cubic(point)[0]
        
    elif method == 'spline_alpha':
        # Hybrid Method: Cubic for (M, H), Cubic Spline for Alpha
        # 1. Extract the alpha slice for the specific M and H using cubic interpolation
        alpha_axis = axes[2]
        slice_points = np.array([[point[0], point[1], a] for a in alpha_axis])
        alpha_slice_values = interp_cubic(slice_points) # Shape: (len(alpha), len(coeffs))
        
        # 2. Fit a cubic smoothing spline along the alpha axis for each coefficient
        coeffs = np.zeros(len(COEFF_COLUMNS))
        for i in range(len(COEFF_COLUMNS)):
            spline = CubicSpline(alpha_axis, alpha_slice_values[:, i])
            coeffs[i] = spline(point[2])
            
    else:
        raise ValueError("Method must be 'linear', 'cubic', or 'spline_alpha'")

    return AeroCoeffs(*coeffs.tolist())