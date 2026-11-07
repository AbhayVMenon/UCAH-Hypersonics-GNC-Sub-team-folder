import os
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
from scipy.optimize import brentq

# ==============================================================================
# CONSTANTS & INPUTS
# ==============================================================================

TARGET_G_LOAD           = 15.0        # Desired normal load factor to pull [g]                      CHANGE TO DESIRED G-LOAD
MACH_NUMBER             = 6        # Flight Mach number [-]                                         CHANGE TO DESIRED MACH
ALTITUDE_M              = 20000.0    # Flight altitude [m] (30 km max)                              CHANGE TO DESIRED ALTITUDE

MASS_KG                 = 120.0      # Vehicle mass [kg]
L_TOTAL_M               = 2.00       # Total vehicle length [m]
C_REF_M                 = L_TOTAL_M  # Reference length used to nondimensionalize Cm [m]
S_REF_M2                = 0.462243   # Aerodynamic reference area [m^2]
X_CG_M                  = 1.12       # Center of gravity station aft of nose tip [m] (58.5% L)
X_HINGE_M               = 1.85       # Control surface hinge line aft of nose tip [m]

DEFAULT_DCM_DDELTA_NOSE = -0.015    # Nose-referenced dCm/ddelta [per deg]

ALPHA_MIN_DEG           = -20.0      # Hard minimum angle of attack limit [deg]
ALPHA_MAX_DEG           = 20.0       # Hard maximum angle of attack limit [deg]
DELTA_LIMIT_DEG         = 25.0       # Physical actuator deflection limit [deg]
INCLUDE_CN_DELTA        = True       # True: infer CN_delta from dCm_ddelta_nose and X_HINGE_M
AUTO_FIX_CM_NOSE_SIGN   = True       # True: ensures nose-down moment at +CN is negative (-CN * xcp/L)

G_ACCEL                 = 9.80665    # Standard gravity [m/s^2]
GAMMA_AIR               = 1.4        # Ratio of specific heats for air [-]
R_GAS_AIR               = 287.05     # Specific gas constant for dry air [J/(kg*K)]

AERO_DATABASE_FILE      = "harv_aero_database.csv"  # Supports .csv, .xlsx, .xls, or .txt/.dat      CHANGE TO YOUR LOCAL DATABASE FILE

COL_ALIASES = {
    "mach":       ["Mach", "mach", "M", "Mach_num"],
    "alt_m":      ["alt_m", "Altitude_m", "alt", "h_m", "Altitude"],
    "alpha_deg":  ["alpha_deg", "alpha", "AoA_deg", "AoA", "Alpha"],
    "CN":         ["CN", "CN0", "C_N", "CN_base"],
    "CA":         ["CA", "CA0", "C_A", "CA_base"],
    "Cm_nose":    ["Cm_nose", "CM_nose", "Cm0_nose", "Cm", "C_m"],
    "xcp_over_L": ["xcp_over_L", "xcp_L", "XCP_OVER_L"],
    "dCm_ddelta": ["dCm_ddelta_nose", "dCm_ddelta", "Cm_delta", "dCmdelta", "C_m_delta"],
}


# ==============================================================================
# LOCAL ATMOSPHERIC CONDITIONS
# ==============================================================================
def get_atmosphere(h_m):
    """Returns density rho [kg/m^3], pressure P [Pa], temp T [K], and sound speed a [m/s]."""
    if h_m < 0.0 or h_m > 47000.0:
        raise ValueError(f"Altitude {h_m} m is outside valid atmosphere model bounds [0, 47000] m.")

    h = float(h_m)
    if h <= 11000.0:
        T = 288.15 - 0.0065 * h
        P = 101325.0 * (T / 288.15) ** (-G_ACCEL / (-0.0065 * R_GAS_AIR))
    elif h <= 20000.0:
        T = 216.65
        P = 22632.06 * np.exp(-G_ACCEL * (h - 11000.0) / (R_GAS_AIR * 216.65))
    elif h <= 32000.0:
        T = 216.65 + 0.0010 * (h - 20000.0)
        P = 5474.889 * (T / 216.65) ** (-G_ACCEL / (0.0010 * R_GAS_AIR))
    else:
        T = 228.65 + 0.0028 * (h - 32000.0)
        P = 868.0187 * (T / 228.65) ** (-G_ACCEL / (0.0028 * R_GAS_AIR))

    rho = P / (R_GAS_AIR * T)
    a   = float(np.sqrt(GAMMA_AIR * R_GAS_AIR * T))
    return rho, P, T, a


# ==============================================================================
# DATABASE PARSER & INTERPOLATOR
# ==============================================================================
class FlexibleAeroDatabase:
    def __init__(self, df):
        self.df = df.copy()
        self.col_map = self._resolve_columns(self.df)

        M_col   = self.col_map["mach"]
        a_col   = self.col_map["alpha_deg"]
        CN_col  = self.col_map["CN"]
        CA_col  = self.col_map["CA"]
        Cm_col  = self.col_map["Cm_nose"]
        alt_col = self.col_map.get("alt_m", None)
        dCm_col = self.col_map.get("dCm_ddelta", None)

        # Ensure Standard Sign Convention on Cm_nose:
        # Positive normal force aft of the nose tip must produce a nose-down (-Cm_nose) moment.
        if AUTO_FIX_CM_NOSE_SIGN:
            pos_mask = self.df[CN_col] > 0.05
            if pos_mask.any() and (self.df.loc[pos_mask, Cm_col] > 0).mean() > 0.5:
                self.df[Cm_col] = -self.df[Cm_col]

        if dCm_col is None:
            dCm_col = "_dCm_ddelta_default"
            self.df[dCm_col] = DEFAULT_DCM_DDELTA_NOSE

        self.use_3d_alt = (alt_col is not None) and (self.df[alt_col].nunique() > 1)
        self._build_interpolators(M_col, alt_col, a_col, [CN_col, CA_col, Cm_col, dCm_col])

        self.alpha_min_db = float(self.df[a_col].min())
        self.alpha_max_db = float(self.df[a_col].max())

    @staticmethod
    def _resolve_columns(df):
        stripped_cols = {c.strip().lower(): c for c in df.columns}
        resolved = {}
        for key, candidates in COL_ALIASES.items():
            for cand in candidates:
                if cand.strip().lower() in stripped_cols:
                    resolved[key] = stripped_cols[cand.strip().lower()]
                    break
        for req in ["mach", "alpha_deg", "CN", "CA", "Cm_nose"]:
            if req not in resolved:
                raise KeyError(f"Required column '{req}' not found. Available columns: {list(df.columns)}")
        return resolved

    def _build_interpolators(self, M_col, alt_col, a_col, value_cols):
        idx_cols = [M_col, alt_col, a_col] if self.use_3d_alt else [M_col, a_col]
        df_sorted = self.df.drop_duplicates(subset=idx_cols).sort_values(by=idx_cols)

        axes = [np.sort(df_sorted[c].unique().astype(float)) for c in idx_cols]
        expected_size = int(np.prod([len(ax) for ax in axes]))

        self.interps = {}
        if len(df_sorted) == expected_size:
            self.is_structured = True
            shape = tuple(len(ax) for ax in axes)
            for name, col in zip(["CN", "CA", "Cm_nose", "dCm_ddelta"], value_cols):
                grid_vals = df_sorted[col].to_numpy(dtype=float).reshape(shape)
                # bounds_error=True strictly forbids extrapolation outside database limits
                self.interps[name] = RegularGridInterpolator(
                    tuple(axes), grid_vals, method="linear", bounds_error=True
                )
        else:
            self.is_structured = False
            pts = df_sorted[idx_cols].to_numpy(dtype=float)
            for name, col in zip(["CN", "CA", "Cm_nose", "dCm_ddelta"], value_cols):
                self.interps[name] = LinearNDInterpolator(pts, df_sorted[col].to_numpy(dtype=float))

    @classmethod
    def from_file(cls, filepath):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Aerodynamic database file not found: '{filepath}'")

        ext = os.path.splitext(filepath)[1].lower()
        if ext in [".xlsx", ".xls"]:
            df = pd.read_excel(filepath)
        elif ext in [".txt", ".dat"]:
            df = pd.read_csv(filepath, sep=r"\s+")
        else:
            df = pd.read_csv(filepath)
        return cls(df)

    def evaluate_state_at_alpha(self, alpha_deg, mach, alt_m=0.0):
        """
        Queries (Mach, [alt_m], alpha_deg), shifts moments from the nose tip (x=0)
        to X_CG_M = 1.17 m, and solves for the required trim control deflection delta_deg.
        """
        pt = np.array([mach, alt_m, alpha_deg], dtype=float) if self.use_3d_alt else np.array([mach, alpha_deg], dtype=float)

        CN0      = self.interps["CN"](pt).item()
        CA0      = self.interps["CA"](pt).item()
        Cm0_nose = self.interps["Cm_nose"](pt).item()
        dCm_nose = self.interps["dCm_ddelta"](pt).item()

        if any(np.isnan(v) for v in (CN0, CA0, Cm0_nose, dCm_nose)):
            raise ValueError(f"Query point (Mach={mach}, alt_m={alt_m}, alpha_deg={alpha_deg}) is outside the database convex hull.")

        # Nondimensional CG shift from nose tip (1.17 m / 2.00 m = 0.585)
        x_cg_bar = X_CG_M / C_REF_M

        # Normal force effectiveness CN_delta [per deg] inferred from nose moment derivative
        if INCLUDE_CN_DELTA and X_HINGE_M > 0.0:
            CN_delta = -(C_REF_M / X_HINGE_M) * dCm_nose
        else:
            CN_delta = 0.0

        # Shift baseline pitching moment and control derivative from nose tip to CG
        Cm0_cg        = Cm0_nose + CN0 * x_cg_bar
        dCm_ddelta_cg = dCm_nose + CN_delta * x_cg_bar

        # Solve CG pitching moment balance: Cm0_cg + dCm_ddelta_cg * delta = 0
        delta_deg = -Cm0_cg / dCm_ddelta_cg

        # Total Normal & Axial coefficients
        CN_tot = CN0 + CN_delta * delta_deg
        CA_tot = CA0 + abs(CN_delta * delta_deg * float(np.sin(np.deg2rad(delta_deg))))

        # Wind-axis Lift & Drag coefficients
        a_rad  = float(np.deg2rad(alpha_deg))
        CL_tot = CN_tot * float(np.cos(a_rad)) - CA_tot * float(np.sin(a_rad))
        CD_tot = CN_tot * float(np.sin(a_rad)) + CA_tot * float(np.cos(a_rad))

        return {
            "alpha_deg":     float(alpha_deg),
            "delta_deg":     float(delta_deg),
            "CL_tot":        float(CL_tot),
            "CD_tot":        float(CD_tot),
            "CN0":           float(CN0),
            "CA0":           float(CA0),
            "CN_tot":        float(CN_tot),
            "CA_tot":        float(CA_tot),
            "Cm0_nose":      float(Cm0_nose),
            "Cm0_cg":        float(Cm0_cg),
            "xcp_over_L":    float(-Cm0_nose / CN0) if abs(CN0) > 1e-6 else float("nan"),
            "dCm_ddelta_cg": float(dCm_ddelta_cg),
        }


# ==============================================================================
# SOLVER
# ==============================================================================
def solve_single_point_trim(n_load_g, mach, alt_m, aero_db):
    rho, P_inf, T_inf, a_sound = get_atmosphere(alt_m)
    V_ms     = mach * a_sound
    q_inf_Pa = 0.5 * rho * (V_ms ** 2)
    W_N      = MASS_KG * G_ACCEL

    def lift_residual(a_deg):
        state = aero_db.evaluate_state_at_alpha(a_deg, mach, alt_m)
        return (q_inf_Pa * S_REF_M2 * state["CL_tot"]) - (n_load_g * W_N)

    # Strictly enforce the intersection of user alpha limits and database alpha bounds
    a_low  = max(ALPHA_MIN_DEG, aero_db.alpha_min_db)
    a_high = min(ALPHA_MAX_DEG, aero_db.alpha_max_db)

    res_low  = lift_residual(a_low)
    res_high = lift_residual(a_high)

    if res_low * res_high > 0:
        n_min_achievable = (res_low + n_load_g * W_N) / W_N
        n_max_achievable = (res_high + n_load_g * W_N) / W_N
        raise ValueError(
            f"Target load factor {n_load_g:.2f} g cannot be achieved within valid alpha bounds "
            f"[{a_low:.1f}, {a_high:.1f}] deg at Mach {mach}, alt {alt_m} m. "
            f"Achievable load factor range is [{n_min_achievable:.2f} g, {n_max_achievable:.2f} g]."
        )

    alpha_trim_deg = brentq(lift_residual, a_low, a_high, xtol=1e-7)

    trimmed = aero_db.evaluate_state_at_alpha(alpha_trim_deg, mach, alt_m)
    trimmed.update({
        "n_load_g":  n_load_g,
        "mach":      mach,
        "alt_m":     alt_m,
        "vel_ms":    V_ms,
        "q_inf_kPa": q_inf_Pa / 1000.0,
        "CL_req":    (n_load_g * W_N) / (q_inf_Pa * S_REF_M2),
        "L_over_D":  trimmed["CL_tot"] / trimmed["CD_tot"] if trimmed["CD_tot"] != 0 else 0.0,
        "saturated": abs(trimmed["delta_deg"]) > DELTA_LIMIT_DEG,
    })
    return trimmed


# ==============================================================================
# FILE EXECUTION
# ==============================================================================
if __name__ == "__main__":
    aero_db = FlexibleAeroDatabase.from_file(AERO_DATABASE_FILE)

    res = solve_single_point_trim(
        n_load_g = TARGET_G_LOAD,
        mach     = MACH_NUMBER,
        alt_m    = ALTITUDE_M,
        aero_db  = aero_db
    )

    print("=" * 64)
    print(f"UNPOWERED HYPERSONIC TRIM SOLUTION (L = {L_TOTAL_M} m, x_cg = {X_CG_M} m)")
    print("=" * 64)
    print(f"  Database Mode        : {'3D (Mach, alt_m, alpha_deg)' if aero_db.use_3d_alt else '2D (Mach, alpha_deg)'}")
    print(f"  Target Maneuver Load : {res['n_load_g']:8.2f} g")
    print(f"  Flight Condition     : Mach {res['mach']:.2f} at {res['alt_m']/1000:.1f} km ({res['vel_ms']:.1f} m/s)")
    print(f"  Dynamic Pressure     : {res['q_inf_kPa']:8.2f} kPa")
    print("-" * 64)
    print(f"  REQUIRED DEFLECTION  : {res['delta_deg']:+8.3f} deg  (Actuator Limit: +/-{DELTA_LIMIT_DEG} deg)")
    print(f"  Trim Angle of Attack : {res['alpha_deg']:+8.3f} deg")
    print("-" * 64)
    print(f"  Required Lift Coeff  : {res['CL_req']:8.4f}  (Trimmed L/D = {res['L_over_D']:.2f})")
    print(f"  Baseline CN0 / CA0   : {res['CN0']:8.4f} / {res['CA0']:.4f}")
    print(f"  Total CN / Total CA  : {res['CN_tot']:8.4f} / {res['CA_tot']:.4f}")
    print(f"  x_cp / L (from table): {res['xcp_over_L']:8.4f}  (x_cg / L = {X_CG_M/L_TOTAL_M:.4f})")
    print(f"  Cm0 (Nose -> CG)     : {res['Cm0_nose']:+8.4f} -> {res['Cm0_cg']:+.4f}")
    print(f"  dCm/ddelta at CG     : {res['dCm_ddelta_cg']:+8.5f} /deg")
    print("=" * 64)