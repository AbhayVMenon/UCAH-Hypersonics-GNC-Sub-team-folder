"""Static trim with a linear CG-moment derivative for ONE bottom flap.

Bottom-pair incremental Cm = 2 * CM_DELTA_BOTTOM_EACH_PER_DEG * delta_bottom_deg.
Linearity and matching operating conditions are assumed across the search domain.
The signed derivative is per degree, uses the vehicle's coefficient references,
and is referenced to zero deployment (closed surfaces are in the baseline).
Top moments, both pairs' forces, and hinge torques retain the pressure model.
The derivative does not independently determine lift, drag, or hinge torque.

Provide --cm-delta-bottom-per-deg YOUR_SIGNED_VALUE on every example below,
or set CM_DELTA_BOTTOM_EACH_PER_DEG in this file. No value is guessed.

Pitch-only quasi-steady trim estimate; CG-referenced aerodynamic database.

Requires numpy, pandas, scipy (openpyxl for xlsx). Run:
    python defAngle3.py --database shadow_comparison.csv
    python defAngle3.py --database shadow_comparison.csv --flight-path-angle 0

Requirement: 15 g TURN = 15*g acceleration perpendicular to velocity.
Lift load factor and total kinematic magnitude are reported separately, not constrained.
Either pitch direction is allowed by default. Supply --direction to restrict it.
Omitting --flight-path-angle permits any flight-path orientation; selected angle is reported.
All four flaps are available: mirrored left/right members move together,
with independent top and bottom deployment. Flap/body interference between
deployed pairs is not predicted; their modeled increments are added.
Joint trim minimizes the largest
pair deployment. Use --allocation one_pair for the earlier restricted model.
Freestream Mach is restricted to 3-8 (local supersonic flow can be slower).
Nominal uniform pressure runs by default; --sensitivity opts in to variations.
Use --maximize --max-deflection 45 to maximize load at a fixed flight condition.
Use --capability to compare loads at advisory/model limits.
Use --envelope --g 15 --max-deflection 45 --envelope-csv turn_envelope.csv
for a nominal-pressure grid study over the database Mach/altitude coverage.
Optional --max-q-kpa and --max-hinge-torque constrain that study.
No stability analysis, pitch-rate aerodynamics, or trajectory propagation.
Local flow and pressure-distribution cases are assumptions until calibrated.
Sensitivity ranges are not statistical confidence intervals or accuracy claims.
Baseline coefficients are used unchanged. Changing S_REF_M2 models body-area
scaling at fixed coefficients; changing C_REF_M scales baseline moment length.
Database coefficient conventions and CG reference must match these inputs.
"""
import os
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
from scipy.optimize import brentq, least_squares, minimize_scalar, minimize
from functools import lru_cache
from dataclasses import dataclass, replace
from pathlib import Path
import argparse

# ==============================================================================
# 1. CONSTANTS & LIFTING-BODY GEOMETRY INPUTS
# ==============================================================================

# --- Target Maneuver Inputs ---
TARGET_G_MAGNITUDE          = 15.0              # Required normal kinematic acceleration |a_perp|/g for turning
MANEUVER_DIRECTION          = "either"          # either, nose_up (+normal), nose_down (-normal)
MACH_NUMBER                 = 5.0               # Flight Mach number [-]
ALTITUDE_M                  = 12000.0           # Flight altitude [m]

# --- Lifting Body Mass & Reference Geometry ---
MASS_KG                     = 101.1849             # Vehicle mass [kg]
L_TOTAL_M                   = 2.00              # Total vehicle length [m]
C_FLAP_M                    = 0.280             # Authoritative flap chord [m], not mm
X_HINGE_M                   = L_TOTAL_M - C_FLAP_M  # Flap ends at vehicle trailing edge

C_REF_M                     = 1.57921038497254  # Reference length used to nondimensionalize Cm [m]
S_REF_M2                    = 0.462243          # Body aerodynamic area scale [m^2], fixed-coefficient approximation
X_CG_M                      = 1.09              # Center of gravity station aft of nose tip [m]
# Use table coefficients directly: N=q*S_REF*CN and M=q*S_REF*C_REF*Cm.
# Area/length studies hold coefficients fixed as a geometry-scaling approximation.
# Flap dimensions are independent; their dimensional loads do not scale with S_REF.
# Table moment reference must match X_CG_M; no normalization or CG conversion.

# --- Outward-Only OML Body Flap Geometry (0 deg = Flush with OML, +deg = Open Outward) ---
# phi measured from HORIZONTAL [deg] (0 = flat upper/lower deck, >0 = canted sidewall/chines)
S_FLAP_TOP_EACH_M2          = C_FLAP_M * 0.0625     # Planform area of ONE top body flap [m^2]     (UPDATE FROM CAD)
PHI_TOP_DEG                 = 44.7                # Top OML dihedral cant angle [deg]            (UPDATE FROM CAD)
Z_HINGE_TOP_M               = 0.05 - 0.014765              # Top hinge vertical offset ABOVE CG [m]       (UPDATE FROM CAD)
ETA_Q_TOP                   = 1.00              # Local dynamic pressure ratio q_top / q_inf [-]

S_FLAP_BOT_EACH_M2          = C_FLAP_M * 0.0735     # Planform area of ONE bottom body flap [m^2]  (UPDATE FROM CAD)
PHI_BOT_DEG                 = 17.9                # Bottom OML anhedral cant angle [deg]         (UPDATE FROM CAD)
Z_HINGE_BOT_M               = -(0.014765 + 0.017)            # Bottom hinge vertical offset BELOW CG [m]    (UPDATE FROM CAD)
ETA_Q_BOT                   = 1.00              # Local dynamic pressure ratio q_bot / q_inf [-]

CP_CAVITY                   = 0.0               # Inner sheltered flap cavity Cp (P_cav = P_inf -> 0.0)
DELTA_BIAS_OUT_DEG          = 0.0               # Minimum outward pre-deployment on all flaps [deg] (0 = flush)
CSV_INCLUDES_CLOSED_OML     = True              # True: subtracts delta_out=0 flush OML pressure already in CSV

# --- Hardware Actuator Limit ---
DELTA_ACTUATOR_MAX_DEG      = 45.0              # Max outward hardware hinge travel [deg] (0 to +DELTA_MAX)

G_ACCEL                     = 9.80665
GAMMA_AIR                   = 1.4
R_GAS_AIR                   = 287.05

AERO_DATABASE_FILE          = "shadow_comparison.csv"

# Hard aerodynamic domain; further restricted to database support. No extrapolation.
ALPHA_MIN_DEG, ALPHA_MAX_DEG = -16.0, 16.0
FLIGHT_PATH_ANGLE_DEG = None  # None: choose a favorable flight-path orientation and report it.
# Any turn profile is allowed. Specifying --flight-path-angle fixes that orientation.
# Pitch-only normal acceleration = L/m - g*cos(gamma); drag is not turning acceleration.
DELTA_MODEL_MAX_DEG = 89.0    # Search domain, independent of advisory hardware travel.
RUN_SENSITIVITY = False  # Uniform nominal run by default; opt in with --sensitivity.
MACH_MIN, MACH_MAX = 3.0, 8.0
# Net pressure shape dp(x)=mean_dp*(1+k*(2*x/c-1)); zero means uniform.
PRESSURE_GRADIENT_TOP = 0.0
PRESSURE_GRADIENT_BOT = 0.0
CONTROL_ALLOCATION = "both_pairs"  # both_pairs (four flaps) or one_pair comparison
# Pitch-only symmetry: left/right top share delta_top; left/right bottom share delta_bot.
# Joint allocation minimizes max(delta_top, delta_bot), with advisory hardware travel.
# Local upstream states are INPUTS, not inferred from the baseline force table.
# None Mach means freestream Mach; pressure ratio is p_edge / p_infinity.
LOCAL_FLOW = {
    "top": {"mach": None, "pressure_ratio": 1.0, "alpha_offset_deg": 0.0},
    "bot": {"mach": None, "pressure_ratio": 1.0, "alpha_offset_deg": 0.0},
}
# ETA_Q_TOP/BOT retained above for legacy configuration, but must stay 1.
# Dynamic-pressure ratios now follow q_edge/q_inf = (p_edge/p_inf)*(M_edge/M_inf)**2.


COL_ALIASES = {
    "mach":      ["Mach", "mach", "M", "Mach_num"],
    "alt_m":     ["alt_m", "Altitude_m", "alt", "h_m", "Altitude"],
    "alpha_deg": ["alpha_deg", "alpha", "AoA_deg", "AoA", "Alpha"],
    "beta_deg":  ["beta_deg", "beta", "Sideslip_deg"],
    "CN":        ["CN_on", "CN", "CN0", "C_N", "CN_base"],
    "CA":        ["CA_on", "CA", "CA0", "C_A", "CA_base"],
    "Cm_cg":     ["Cm_on", "Cm_cg", "Cm", "C_m"],
}


# New input for defAngle3; all inherited constants above remain unchanged.
# Signed change in VEHICLE Cm per DEGREE for ONE bottom flap, about CG.
# Set your measured value here, or supply --cm-delta-bottom-per-deg at runtime.
CM_DELTA_BOTTOM_EACH_PER_DEG = 0.0004359


# ==============================================================================
# 2. ATMOSPHERE & COMPRESSIBLE OML PANEL AERODYNAMICS
# ==============================================================================
def get_atmosphere(h_m):
    if h_m < 0.0 or h_m > 47000.0:
        raise ValueError(f"Altitude {h_m} m outside valid bounds [0, 47000] m.")
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
    a = float(np.sqrt(GAMMA_AIR * R_GAS_AIR * T))
    return rho, P, T, a


def rayleigh_pitot_cp_max(mach, gamma=GAMMA_AIR):
    M2 = mach ** 2
    t1 = (((gamma + 1.0) ** 2 * M2) / (4.0 * gamma * M2 - 2.0 * (gamma - 1.0))) ** (gamma / (gamma - 1.0))
    t2 = (1.0 - gamma + 2.0 * gamma * M2) / (gamma + 1.0)
    return float((2.0 / (gamma * M2)) * (t1 * t2 - 1.0))


def prandtl_meyer_nu(mach, gamma=GAMMA_AIR):
    if mach <= 1.0:
        return 0.0
    g_ratio = np.sqrt((gamma + 1.0) / (gamma - 1.0))
    m_term  = np.sqrt(mach ** 2 - 1.0)
    return float(g_ratio * np.arctan(m_term / g_ratio) - np.arctan(m_term))


@dataclass(frozen=True)
class PressureCase:
    """Coherent pressure assumptions applied to both closed and deployed states.

    cp_scale changes the exterior pressure departure from local static pressure.
    k specifies net-pressure shape: dp(x)=mean_dp*(1+k*(2*x/c-1)).
    Uniform k=0; k=-1/+1 gives CP at c/3 or 2c/3. These are sensitivity
    assumptions, not measured bounds. cavity_cp is referenced to freestream.
    """
    name: str = "nominal"
    cp_scale_top: float = 1.0
    cp_scale_bot: float = 1.0
    k_top: float = 0.0
    k_bot: float = 0.0
    cavity_cp: float = CP_CAVITY
    compression: str = "shock_expansion"  # or modified_newtonian comparison
    expansion: str = "prandtl_meyer"      # or ambient comparison
    flow_angle_offset_deg: float = 0.0


NOMINAL_CASE = PressureCase(k_top=PRESSURE_GRADIENT_TOP, k_bot=PRESSURE_GRADIENT_BOT)
# At the defaults, pressure is uniform on each face and CP is at half chord.
# Different pairs/states can have different pressure MAGNITUDES.
# Deliberately illustrative variations. Replace amplitudes with calibrated bounds.
PRESSURE_CASES = [
    NOMINAL_CASE,
    replace(NOMINAL_CASE, name="exterior_cp_0.75", cp_scale_top=.75, cp_scale_bot=.75),
    replace(NOMINAL_CASE, name="exterior_cp_1.25", cp_scale_top=1.25, cp_scale_bot=1.25),
    replace(NOMINAL_CASE, name="cp_forward", k_top=-1.0, k_bot=-1.0),
    replace(NOMINAL_CASE, name="cp_aft", k_top=1.0, k_bot=1.0),
    replace(NOMINAL_CASE, name="modified_newtonian", compression="modified_newtonian"),
    replace(NOMINAL_CASE, name="no_expansion_suction", expansion="ambient"),
    replace(NOMINAL_CASE, name="cavity_cp_plus_0.02", cavity_cp=.02),
    replace(NOMINAL_CASE, name="cavity_cp_minus_0.02", cavity_cp=-.02),
    replace(NOMINAL_CASE, name="local_flow_minus_2deg", flow_angle_offset_deg=-2.),
    replace(NOMINAL_CASE, name="local_flow_plus_2deg", flow_angle_offset_deg=2.),
    replace(NOMINAL_CASE, name="top_weak_bottom_strong", cp_scale_top=.75, cp_scale_bot=1.25),
    replace(NOMINAL_CASE, name="top_strong_bottom_weak", cp_scale_top=1.25, cp_scale_bot=.75),
]


def validate_freestream_mach(mach):
    if not np.isfinite(mach) or not MACH_MIN <= mach <= MACH_MAX:
        raise ValueError(f"Freestream Mach must be within [{MACH_MIN:g}, {MACH_MAX:g}]; got {mach}.")


@lru_cache(maxsize=128)
def attached_shock_limit(mach, gamma=GAMMA_AIR):
    """Shock angle and maximum attached turning angle for perfect-gas 2D flow."""
    mu = np.arcsin(1.0 / mach)
    def turning(beta):
        return np.arctan(2 / np.tan(beta) * (mach**2 * np.sin(beta)**2 - 1)
                         / (mach**2 * (gamma + np.cos(2*beta)) + 2))
    fit = minimize_scalar(lambda b: -turning(b), bounds=(mu, np.pi/2), method="bounded")
    return float(fit.x), float(turning(fit.x))


def pressure_cp(theta_rad, mach, case=NOMINAL_CASE, gamma=GAMMA_AIR):
    """Return local-reference Cp and a model-regime label.

    Attached compression: exact weak oblique shock. Detached compression:
    explicitly flagged modified-Newtonian fallback, not a separation solution.
    Expansion: ideal Prandtl-Meyer, without the old arbitrary suction fade.
    """
    if mach <= 1.0:
        raise ValueError("Pressure model requires supersonic local Mach.")
    if abs(theta_rad) < 1e-12:
        return 0.0, "parallel"
    if theta_rad > 0:
        if case.compression == "modified_newtonian":
            return rayleigh_pitot_cp_max(mach, gamma)*np.sin(theta_rad)**2, "modified_newtonian"
        if case.compression != "shock_expansion":
            raise ValueError("Unknown compression model.")
        beta_peak, theta_max = attached_shock_limit(mach, gamma)
        if theta_rad > theta_max:
            return rayleigh_pitot_cp_max(mach, gamma)*np.sin(theta_rad)**2, "detached_fallback"
        def residual(beta):
            return (2/np.tan(beta)*(mach**2*np.sin(beta)**2-1)
                    /(mach**2*(gamma+np.cos(2*beta))+2) - np.tan(theta_rad))
        beta = brentq(residual, np.arcsin(1/mach), beta_peak, xtol=1e-12)
        return 4/(gamma+1)*(np.sin(beta)**2-1/mach**2), "attached_shock"
    if case.expansion == "ambient":
        return 0.0, "no_expansion_suction"
    if case.expansion != "prandtl_meyer":
        raise ValueError("Unknown expansion model.")
    nu_target = prandtl_meyer_nu(mach, gamma) - theta_rad
    nu_max = .5*np.pi*(np.sqrt((gamma+1)/(gamma-1))-1)
    if nu_target >= nu_max:
        return -2/(gamma*mach**2), "vacuum_limit"
    # Dynamically bracket the inverse rather than clipping at Mach 1000.
    upper = max(2*mach, 2.)
    while prandtl_meyer_nu(upper, gamma) < nu_target:
        upper *= 2
    m2 = brentq(lambda m: prandtl_meyer_nu(m, gamma)-nu_target, mach, upper)
    p_ratio = ((1+.5*(gamma-1)*mach**2)/(1+.5*(gamma-1)*m2**2))**(gamma/(gamma-1))
    return 2/(gamma*mach**2)*(p_ratio-1), "attached_expansion"


def surface_face_cp(theta_rad, mach, cp_max=None, gamma=GAMMA_AIR):
    """Compatibility wrapper; cp_max is recomputed from local Mach."""
    return pressure_cp(theta_rad, mach, NOMINAL_CASE, gamma)[0]


def flap_pair_raw(pair_name, alpha_deg, delta_deg, mach, case=NOMINAL_CASE):
    """Integrate a rectangular flat pair's linear net-pressure distribution.

    Coordinates: x aft, z up; positive pitching moment nose up.
    Analytical integration is exact for the ASSUMED linear distribution.
    It does not predict separation, shock impingement, or body interference.
    """
    if pair_name not in ("top", "bot"):
        raise ValueError("pair_name must be top or bot.")
    if not 0 <= delta_deg <= DELTA_MODEL_MAX_DEG:
        raise ValueError("Deflection outside aerodynamic search domain.")
    validate_freestream_mach(mach)
    top = pair_name == "top"
    area = S_FLAP_TOP_EACH_M2 if top else S_FLAP_BOT_EACH_M2
    phi = np.deg2rad(PHI_TOP_DEG if top else PHI_BOT_DEG)
    z_h = Z_HINGE_TOP_M if top else Z_HINGE_BOT_M
    k = case.k_top if top else case.k_bot
    scale = case.cp_scale_top if top else case.cp_scale_bot
    if abs(k) > 1 or scale < 0:
        raise ValueError("Require |k| <= 1 and nonnegative exterior Cp scale.")
    flow = LOCAL_FLOW[pair_name]
    me = mach if flow["mach"] is None else float(flow["mach"])
    pe_ratio = float(flow["pressure_ratio"])
    if pe_ratio <= 0:
        raise ValueError("Local static pressure ratio must be positive.")
    a = np.deg2rad(alpha_deg + flow["alpha_offset_deg"] + case.flow_angle_offset_deg)
    d = np.deg2rad(delta_deg)
    sign = -1 if top else 1
    sin_theta = np.cos(a)*np.sin(d) + sign*np.sin(a)*np.cos(d)*np.cos(phi)
    theta = np.arcsin(np.clip(sin_theta, -1, 1))
    cp_local, regime = pressure_cp(theta, me, case)
    q_ratio = pe_ratio*(me/mach)**2
    # p_outer = p_edge + q_edge*Cp_local; inner uses global Cp.
    cp_outer = (pe_ratio-1)*2/(GAMMA_AIR*mach**2) + q_ratio*scale*cp_local
    cp_net = cp_outer-case.cavity_cp
    s_cp = C_FLAP_M*(.5+k/6)
    c_panel = 2*area/S_REF_M2*cp_net
    cn = sign*c_panel*np.cos(d)*np.cos(phi)
    ca = c_panel*np.sin(d)
    cm = (z_h*ca-(X_HINGE_M-X_CG_M)*cn-sign*s_cp*np.cos(phi)*c_panel)/C_REF_M
    return {"CN": cn, "CA": ca, "Cm": cm, "cp_net": cp_net,
            "cp_local": cp_local, "cp_outer": cp_outer,
            "cp_distance_m": s_cp, "regime": regime,
            "area_each_m2": area, "q_ratio": q_ratio,
            "theta_deg": float(np.rad2deg(theta)),
            "attached_theta_max_deg": float(np.rad2deg(attached_shock_limit(me)[1])),
            "pressure_gradient_k": k}


def bottom_pair_cm_increment(delta_bottom_deg, derivative_each_per_deg=None):
    """Two symmetric flaps: Delta Cm = 2 * (dCm/dDelta)_each * Delta_deg.

    The baseline supplies the zero-deployment moment. The derivative replaces
    the bottom pair's pressure-derived vehicle moment; it is not added to it.
    Forces and gross per-flap hinge moments still use the pressure model.
    """
    slope=(CM_DELTA_BOTTOM_EACH_PER_DEG if derivative_each_per_deg is None
           else derivative_each_per_deg)
    if slope is None or not np.isfinite(float(slope)):
        raise ValueError('Supply a finite signed Cm derivative for ONE bottom flap per degree '
                         'using --cm-delta-bottom-per-deg or CM_DELTA_BOTTOM_EACH_PER_DEG.')
    if not CSV_INCLUDES_CLOSED_OML:
        raise ValueError('The linear bottom moment requires a baseline including closed surfaces '
                         '(CSV_INCLUDES_CLOSED_OML=True).')
    return float(2.0*float(slope)*float(delta_bottom_deg))


def evaluate_oml_pair_increment(pair_name, alpha_deg, delta_out_deg, mach, cp_max=None,
                                case=NOMINAL_CASE):
    raw = flap_pair_raw(pair_name, alpha_deg, delta_out_deg, mach, case)
    closed = flap_pair_raw(pair_name, alpha_deg, 0., mach, case) if CSV_INCLUDES_CLOSED_OML else {}
    increments=tuple(raw[name]-closed.get(name, 0.) for name in ("CN", "CA", "Cm"))
    if pair_name == "bot":
        return increments[0], increments[1], bottom_pair_cm_increment(delta_out_deg)
    return increments


class FlexibleAeroDatabase:
    def __init__(self, df):
        self.df = df.copy()
        self.col_map = self._resolve_columns(self.df)
        if "beta_deg" in self.col_map:
            b_col = self.col_map["beta_deg"]
            self.df = self.df[np.isclose(self.df[b_col], 0.0, atol=1e-5)].copy()

        M_col   = self.col_map["mach"]
        a_col   = self.col_map["alpha_deg"]
        CN_col  = self.col_map["CN"]
        CA_col  = self.col_map["CA"]
        Cm_col  = self.col_map["Cm_cg"]
        alt_col = self.col_map.get("alt_m", None)

        if self.df.empty:
            raise ValueError("No zero-sideslip rows remain in the database.")

        self.use_3d_alt = (alt_col is not None) and (self.df[alt_col].nunique() > 1)
        self.mach_bounds = (float(self.df[M_col].min()), float(self.df[M_col].max()))
        if self.use_3d_alt:
            self.alt_bounds = (float(self.df[alt_col].min()), float(self.df[alt_col].max()))
        else:
            self.alt_bounds = None

        self.alpha_axis = np.sort(self.df[a_col].unique().astype(float))
        self.alpha_min_db = float(self.alpha_axis[0])
        self.alpha_max_db = float(self.alpha_axis[-1])
        self._build_interpolators(M_col, alt_col, a_col, [CN_col, CA_col, Cm_col])

    @staticmethod
    def _resolve_columns(df):
        # Aerodynamic capitalization is meaningful: CN is normal force,
        # Cn is yawing moment; CL and Cl similarly differ. Prefer exact case.
        exact = {c.strip(): c for c in df.columns}
        resolved = {}
        for key, candidates in COL_ALIASES.items():
            match = next((exact[cand] for cand in candidates if cand in exact), None)
            if match is not None:
                resolved[key] = match
                continue
            # Permit case-insensitive aliases only when there is one unambiguous column.
            for cand in candidates:
                matches = [c for c in df.columns if c.strip().lower() == cand.lower()]
                if len(matches) > 1:
                    raise ValueError(f"Ambiguous columns for {key}: {matches}; use exact alias names.")
                if matches:
                    resolved[key] = matches[0]
                    break
        for req in ["mach", "alpha_deg", "CN", "CA", "Cm_cg"]:
            if req not in resolved:
                raise KeyError(f"Required column '{req}' not found. Available: {list(df.columns)}")
        return resolved

    def _build_interpolators(self, M_col, alt_col, a_col, value_cols):
        idx_cols  = [M_col, alt_col, a_col] if self.use_3d_alt else [M_col, a_col]
        if self.df.duplicated(subset=idx_cols).any():
            raise ValueError("Duplicate database coordinates: resolve configurations before interpolation.")
        df_sorted = self.df.sort_values(by=idx_cols)
        if not np.isfinite(df_sorted[idx_cols + value_cols].to_numpy(dtype=float)).all():
            raise ValueError("Database contains non-finite coordinates or coefficients.")
        axes      = [np.sort(df_sorted[c].unique().astype(float)) for c in idx_cols]
        self.interps = {}
        self.is_structured = (len(df_sorted) == int(np.prod([len(ax) for ax in axes])))

        if self.is_structured:
            shape = tuple(len(ax) for ax in axes)
            for name, col in zip(["CN", "CA", "Cm_cg"], value_cols):
                # Never extrapolate aerodynamic coefficients.
                self.interps[name] = RegularGridInterpolator(
                    tuple(axes), df_sorted[col].to_numpy(dtype=float).reshape(shape),
                    method="linear", bounds_error=True
                )
        else:
            pts = df_sorted[idx_cols].to_numpy(dtype=float)
            for name, col in zip(["CN", "CA", "Cm_cg"], value_cols):
                self.interps[name] = LinearNDInterpolator(pts, df_sorted[col].to_numpy(dtype=float))

    @classmethod
    def from_file(cls, filepath):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Database file not found: '{filepath}'")
        ext = os.path.splitext(filepath)[1].lower()
        df  = pd.read_excel(filepath) if ext in [".xlsx", ".xls"] else (
              pd.read_csv(filepath, sep=r"\s+") if ext in [".txt", ".dat"] else pd.read_csv(filepath))
        return cls(df)

    def get_baseline(self, alpha_deg, mach, alt_m=0.0):
        validate_freestream_mach(mach)
        if not max(ALPHA_MIN_DEG, self.alpha_min_db) <= alpha_deg <= min(ALPHA_MAX_DEG, self.alpha_max_db):
            raise ValueError("AoA outside configured/database bounds; extrapolation is disabled.")
        if not self.mach_bounds[0] <= mach <= self.mach_bounds[1]:
            raise ValueError(f"Mach {mach} outside database bounds {self.mach_bounds}.")
        if self.use_3d_alt and not self.alt_bounds[0] <= alt_m <= self.alt_bounds[1]:
            raise ValueError(f"Altitude {alt_m} outside database bounds {self.alt_bounds}.")
        pt = np.array([mach, alt_m, alpha_deg] if self.use_3d_alt else [mach, alpha_deg])
        values = tuple(float(self.interps[name](pt).item()) for name in ("CN", "CA", "Cm_cg"))
        if not np.isfinite(values).all():
            raise ValueError("Requested condition outside the database interpolation support.")
        # Fixed table coefficients: body-area changes scale dimensional loads.
        return values


def evaluate_state(alpha_deg, delta_top_deg, delta_bot_deg, mach, alt_m, aero_db,
                   case=NOMINAL_CASE):
    cn0, ca0, cm0 = aero_db.get_baseline(alpha_deg, mach, alt_m)
    top = evaluate_oml_pair_increment("top", alpha_deg, delta_top_deg, mach, case=case)
    bot = evaluate_oml_pair_increment("bot", alpha_deg, delta_bot_deg, mach, case=case)
    cn, ca, cm = np.array((cn0, ca0, cm0)) + np.array(top) + np.array(bot)
    a = np.deg2rad(alpha_deg)
    return {"alpha_deg": float(alpha_deg), "delta_top_deg": float(delta_top_deg),
            "delta_bot_deg": float(delta_bot_deg), "CN0": cn0, "CA0": ca0,
            "Cm0_cg": cm0, "CN_tot": cn, "CA_tot": ca, "Cm_total": cm,
            "dCm_top": top[2], "dCm_bot": bot[2], "dCm_flap": top[2]+bot[2],
            "CL_tot": cn*np.cos(a)-ca*np.sin(a),
            "CD_tot": cn*np.sin(a)+ca*np.cos(a)}


def acceleration_metrics(cl, cd, q_inf, g_mag, direction_sign, flight_path_angle_deg):
    """Turning at CG: a_perp/g = L/(mg)-cos(gamma), signed in the pitch plane.

    direction_sign: +1 requests positive turning, -1 negative, 0 either.
    No separate lift or total-acceleration requirement is imposed. Drag is
    reported in tangential acceleration but cannot satisfy turning demand.
    If gamma is unspecified, any pitch-plane turn orientation is allowed:
    choose gamma=180 deg for positive normal acceleration or 0 for negative.
    This is a best-orientation capability, not capability at every orientation.
    Uniform gravity contributes at most 1 g; a 15 g turn still needs aerodynamics.
    """
    n_lift=float(q_inf*S_REF_M2*cl/(MASS_KG*G_ACCEL))
    d_ratio=float(q_inf*S_REF_M2*cd/(MASS_KG*G_ACCEL))
    aero_magnitude_g=float(np.hypot(n_lift,d_ratio))
    lower_bound=abs(aero_magnitude_g-1.)
    upper_bound=aero_magnitude_g+1.
    if direction_sign not in (-1,0,1):
        raise ValueError('direction_sign must be -1, 0 (either) or +1.')
    selected=flight_path_angle_deg is None
    if selected:
        best_sign=direction_sign if direction_sign else (1 if n_lift>=0 else -1)
        gamma_deg=180. if best_sign>0 else 0.
    else:
        gamma_deg=float(flight_path_angle_deg)
    gamma=np.deg2rad(gamma_deg)
    tangential=-d_ratio-np.sin(gamma)
    normal=n_lift-np.cos(gamma)
    kin=float(np.hypot(tangential,normal))
    turning=abs(normal) if direction_sign==0 else direction_sign*normal
    return {"n_load_g":n_lift,"n_kinematic_g":kin,
            "n_normal_g":float(normal),"n_tangential_g":float(tangential),
            "turning_capacity_g":float(turning),"turning_magnitude_g":float(abs(normal)),
            "turning_magnitude_m_s2":float(abs(normal)*G_ACCEL),
            "effective_flight_path_angle_deg":gamma_deg,
            "flight_path_angle_was_selected":selected,
            "n_kinematic_lower_bound_g":lower_bound,"n_kinematic_upper_bound_g":upper_bound,
            "n_kinematic_magnitude_m_s2":kin*G_ACCEL,
            "n_kinematic_range_m_s2":(lower_bound*G_ACCEL,upper_bound*G_ACCEL),
            "requirement_margin_g":float(turning-g_mag)}


class TrimSearchError(ValueError):
    """Finite search failed; diagnostics describe sampled trim states, not a proof."""
    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


def scan_trim_curves(mach, alt_m, aero_db, g_mag, direction_sign,
                     flight_path_angle_deg=None, case=NOMINAL_CASE,
                     alpha_samples=161, delta_samples=181):
    """Bracket Cm=0 in deployment first, then bracket the acceleration boundary.

    Includes rising/falling authority branches, rejects sign changes across pressure
    jumps, and verifies both final residuals. Finite resolution can still miss
    tangent or very narrow solution branches; ranges are sampled, not certified.
    Hardware travel is advisory. One pair remains at bias, as in the main solver.
    """
    lo=max(ALPHA_MIN_DEG,aero_db.alpha_min_db)
    hi=min(ALPHA_MAX_DEG,aero_db.alpha_max_db)
    rho,_,_,sound=get_atmosphere(alt_m)
    q_inf=.5*rho*(mach*sound)**2
    bias=float(DELTA_BIAS_OUT_DEG)
    alphas=np.unique(np.r_[np.linspace(lo,hi,alpha_samples),
                          aero_db.alpha_axis[(aero_db.alpha_axis>=lo)&(aero_db.alpha_axis<=hi)]])
    deltas=np.linspace(bias,DELTA_MODEL_MAX_DEG,delta_samples)
    candidates=[]
    feasible_samples=[]
    sampled=[]
    eval_errors=[]
    for pair in ('top','bot'):
        other='bot' if pair=='top' else 'top'
        cache={}
        def trims(alpha):
            key=float(alpha)
            if key in cache:return cache[key]
            try:
                cn0,ca0,cm0=aero_db.get_baseline(key,mach,alt_m)
            except ValueError as exc:
                eval_errors.append(str(exc));cache[key]=[];return []
            inactive=np.array(evaluate_oml_pair_increment(other,key,bias,mach,case=case))
            closed=flap_pair_raw(pair,key,0.,mach,case) if CSV_INCLUDES_CLOSED_OML else {}
            def moment(d):
                return cm0+inactive[2]+evaluate_oml_pair_increment(pair,key,d,mach,case=case)[2]
            ms=np.array([moment(d) for d in deltas])
            roots=[]
            for i in range(len(deltas)-1):
                if abs(ms[i])<=1e-10:roots.append(float(deltas[i]))
                if ms[i]*ms[i+1]<0:
                    root=brentq(moment,deltas[i],deltas[i+1],xtol=1e-10)
                    # A bracket around a pressure-model jump is not necessarily a root.
                    if abs(moment(root))<=1e-7:roots.append(float(root))
            if abs(ms[-1])<=1e-10:roots.append(float(deltas[-1]))
            states=[]
            for d in sorted(set(round(r,9) for r in roots)):
                dt,db=(d,bias) if pair=='top' else (bias,d)
                state=evaluate_state(key,dt,db,mach,alt_m,aero_db,case)
                state.update(acceleration_metrics(state['CL_tot'],state['CD_tot'],q_inf,
                                                   g_mag,direction_sign,flight_path_angle_deg))
                if state['CD_tot'] < -1e-9:continue
                state['active_pair']=pair.upper() if d>bias+1e-6 else 'NONE'
                state['solver_nfev']=None
                state['search_method']='bracketed_trim_curve'
                states.append(state)
            cache[key]=states
            return states
        previous=[]
        last_alpha=None
        for alpha in alphas:
            current=trims(alpha)
            sampled.extend(current)
            for st in current:
                if abs(st['requirement_margin_g'])<=1e-5:
                    candidates.append(st)
                elif st['requirement_margin_g']>0 and abs(st['Cm_total'])<=1e-7:
                    feasible_samples.append(st)
                if last_alpha is None:continue
                d=st[f'delta_{pair}_deg']
                # Pair nearby sampled roots to avoid obviously different branches.
                for old in previous:
                    old_d=old[f'delta_{pair}_deg']
                    if abs(d-old_d)>15.:continue
                    if old['requirement_margin_g']*st['requirement_margin_g']>=0:continue
                    def choose(a):
                        options=trims(a)
                        if not options:raise ValueError('Trim branch disappears inside bracket.')
                        predicted=old_d+(d-old_d)*(a-last_alpha)/(alpha-last_alpha)
                        return min(options,key=lambda r:abs(r[f'delta_{pair}_deg']-predicted))
                    try:
                        root=brentq(lambda a:choose(a)['requirement_margin_g'],last_alpha,alpha,xtol=1e-9)
                        result=choose(root)
                    except ValueError:
                        continue
                    if abs(result['Cm_total'])<=1e-7 and abs(result['requirement_margin_g'])<=1e-5:
                        candidates.append(result)
            previous=current
            last_alpha=alpha
    diagnostics={'sampled_trim_states':len(sampled),'alpha_bounds_deg':(lo,hi),
                 'allocation':'one pair active, other at bias',
                 'evaluation_errors':list(dict.fromkeys(eval_errors)),
                 'search_resolution':(len(alphas),len(deltas))}
    if sampled:
        best=max(sampled,key=lambda r:r['requirement_margin_g'])
        diagnostics.update({'sampled_lift_range_g':(min(r['n_load_g'] for r in sampled),
                                                   max(r['n_load_g'] for r in sampled)),
                            'best_sampled_requirement_margin_g':best['requirement_margin_g'],
                            'best_sampled_state':best})
    if not candidates and feasible_samples:
        # Requirements are minima; exceeding them is feasible even if the exact
        # boundary lies outside the allowed domain. Prefer the smallest overshoot.
        best=min(feasible_samples,key=lambda r:r['requirement_margin_g']).copy()
        best['search_method']='bracketed_trim_curve_feasible_sample'
        candidates=[best]
    return candidates,diagnostics


def solve_single_pair_trim(g_mag, direction, mach, alt_m, aero_db,
                            flight_path_angle_deg=FLIGHT_PATH_ANGLE_DEG,
                            case=NOMINAL_CASE):
    """Solve two equalities for alpha and ONE active flap pair at a time.

    Both pairs are tried; the inactive pair stays at bias. Multi-start bounded
    least squares handles non-monotone authority without assuming moment sign.
    Bracketed trim-curve fallback handles local optimizer failures.
    Hardware travel never bounds the solve. Verify moment balance and requirements.
    Search is finite: failure means no solution found, not proof of impossibility.
    """
    lo, hi, bias, sgn, vel, q_inf = _joint_context(
        g_mag, direction, mach, alt_m, aero_db, flight_path_angle_deg, case)
    candidates = []
    failures = []
    diagnostics = {}
    for pair in ("top", "bot"):
        def state_at(x):
            dt, db = (x[1], bias) if pair == "top" else (bias, x[1])
            state = evaluate_state(x[0], dt, db, mach, alt_m, aero_db, case)
            state.update(acceleration_metrics(state["CL_tot"], state["CD_tot"], q_inf,
                                              g_mag, sgn, flight_path_angle_deg))
            return state
        def residual(x):
            state = state_at(x)
            return [state["Cm_total"]/.05, state["requirement_margin_g"]/max(g_mag, 1.)]
        alpha_starts = np.unique(np.r_[np.linspace(lo, hi, 5), np.clip(0., lo, hi)])
        for aa in alpha_starts:
            for dd in (bias+.1*(DELTA_MODEL_MAX_DEG-bias),
                       bias+.45*(DELTA_MODEL_MAX_DEG-bias),
                       bias+.85*(DELTA_MODEL_MAX_DEG-bias)):
                try:
                    fit = least_squares(residual, [aa, dd], bounds=([lo, bias], [hi, DELTA_MODEL_MAX_DEG]),
                                        x_scale=[max(hi-lo, 1.), 45.],
                                        ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=180)
                    state = state_at(fit.x)
                except ValueError as exc:
                    failures.append(str(exc))
                    continue
                # Reject roots near discontinuous model switches or inconsistent tables.
                if abs(state["Cm_total"]) > 1e-7 or abs(state["requirement_margin_g"]) > 1e-5:
                    continue
                if state["CD_tot"] < -1e-9:
                    continue
                if any(abs(state["alpha_deg"]-c["alpha_deg"]) < 1e-4 and
                       abs(state["delta_top_deg"]-c["delta_top_deg"]) < 1e-3 and
                       abs(state["delta_bot_deg"]-c["delta_bot_deg"]) < 1e-3 for c in candidates):
                    continue
                state["active_pair"] = pair.upper() if fit.x[1] > bias+1e-6 else "NONE"
                state["solver_nfev"] = fit.nfev
                state["search_method"] = "multistart_least_squares"
                candidates.append(state)
    if not candidates:
        candidates,diagnostics=scan_trim_curves(mach,alt_m,aero_db,g_mag,sgn,
                                               flight_path_angle_deg,case)
    if not candidates:
        if diagnostics.get('sampled_trim_states',0):
            nlo,nhi=diagnostics['sampled_lift_range_g']
            detail=(f" Sampled signed trimmed L/(mg) range: [{nlo:+.3f},{nhi:+.3f}] g."
                    f" Best requirement margin: {diagnostics['best_sampled_requirement_margin_g']:+.3f} g.")
        else:
            detail=" No Cm=0 states found on the sampled trim curves."
        if diagnostics.get('evaluation_errors'):
            detail+=' Database errors: '+diagnostics['evaluation_errors'][0]
        raise TrimSearchError(
            f"No verified trim solution found for {g_mag:g} g in AoA [{lo:g},{hi:g}] "
            f"and flap [{bias:g},{DELTA_MODEL_MAX_DEG:g}] deg."+detail+
            " Hardware travel was not enforced. Finite search, not a proof of infeasibility.",diagnostics)
    # Retain original preference for low-|alpha|, with deterministic tie-breaks.
    result = min(candidates, key=lambda c: (abs(c["alpha_deg"]),
                                            max(c["delta_top_deg"], c["delta_bot_deg"]))).copy()
    warnings = ["Quasi-steady estimate; local flow and pressure distribution need calibration."]
    if flight_path_angle_deg is None:
        warnings.append("Flight-path angle was free: reported acceleration uses the selected favorable orientation, not all orientations.")
    if max(result["delta_top_deg"], result["delta_bot_deg"]) > DELTA_ACTUATOR_MAX_DEG:
        warnings.append(f"Required deployment exceeds advisory {DELTA_ACTUATOR_MAX_DEG:g} deg hardware travel.")
    for pair in ("top", "bot"):
        raw = flap_pair_raw(pair, result["alpha_deg"], result[f"delta_{pair}_deg"], mach, case)
        closed = flap_pair_raw(pair, result["alpha_deg"], 0., mach, case)
        result[f"{pair}_pressure"] = raw
        # Gross hinge load on ONE deployed flap, not the incremental vehicle load.
        result[f"hinge_torque_{pair}_each_Nm"] = abs(q_inf*raw["area_each_m2"]*raw["cp_net"]*raw["cp_distance_m"])
        regimes = {raw["regime"], closed["regime"]} if CSV_INCLUDES_CLOSED_OML else {raw["regime"]}
        if "detached_fallback" in regimes:
            warnings.append(f"{pair}: detached-compression modified-Newtonian fallback; pressure fidelity uncertain.")
        if abs(raw['theta_deg']-raw['attached_theta_max_deg']) < 1.0:
            warnings.append(f"{pair}: incidence is within 1 deg of attached-shock limit; model-switch sensitivity.")
        if "vacuum_limit" in regimes:
            warnings.append(f"{pair}: ideal expansion reached vacuum limit.")
        if result[f"delta_{pair}_deg"] >= DELTA_MODEL_MAX_DEG-.1:
            warnings.append(f"{pair}: solution near aerodynamic search boundary.")
    if flight_path_angle_deg is not None and sgn*result["n_normal_g"] < 0:
        warnings.append("Turning direction opposes the requested direction.")
    result.update({"mach": mach, "alt_m": alt_m, "vel_ms": vel,
                   "q_inf_kPa": q_inf/1000., "target_g": g_mag,
                   "flight_path_angle_deg": result.get("effective_flight_path_angle_deg",flight_path_angle_deg),
                   "requested_flight_path_angle_deg":flight_path_angle_deg,
                   "pressure_case": case.name, "warnings": warnings,
                   "search_diagnostics": diagnostics,
                   "saturated": max(result["delta_top_deg"], result["delta_bot_deg"]) > DELTA_ACTUATOR_MAX_DEG,
                   "alpha_extrapolated": False,
                   "L_over_D": abs(result["CL_tot"])/result["CD_tot"] if result["CD_tot"] > 0 else None,
                   "candidate_count": len(candidates),
                   "candidate_solutions": [{k: c[k] for k in ("alpha_deg", "delta_top_deg", "delta_bot_deg", "Cm_total")}
                                           for c in candidates]})
    return result

def _joint_context(g_mag, direction, mach, alt_m, aero_db, flight_path_angle_deg, case):
    if direction.lower() not in ("nose_up", "nose_down", "either"):
        raise ValueError("direction must be either, nose_up or nose_down.")
    validate_freestream_mach(mach)
    if not np.isfinite([g_mag, alt_m]).all() or g_mag <= 0:
        raise ValueError("Require finite inputs and positive target g.")
    if flight_path_angle_deg is not None and not np.isfinite(flight_path_angle_deg):
        raise ValueError("Flight-path angle must be finite or None.")
    if ETA_Q_TOP != 1.0 or ETA_Q_BOT != 1.0:
        raise ValueError("Set legacy ETA_Q values to 1; configure LOCAL_FLOW instead.")
    bias = float(DELTA_BIAS_OUT_DEG)
    if not 0 <= bias < DELTA_MODEL_MAX_DEG < 90:
        raise ValueError("Require 0 <= bias < model search limit < 90 degrees.")
    lo = max(ALPHA_MIN_DEG, aero_db.alpha_min_db)
    hi = min(ALPHA_MAX_DEG, aero_db.alpha_max_db)
    if lo >= hi:
        raise ValueError("No nonzero AoA interval shared by database and [-16,16] limits.")
    sgn = 0 if direction.lower() == "either" else (1 if direction.lower() == "nose_up" else -1)
    rho, _, _, sound = get_atmosphere(alt_m)
    vel = mach*sound
    q_inf = .5*rho*vel**2
    # Expose configuration errors before entering local optimization attempts.
    if min(S_REF_M2,C_REF_M,MASS_KG,C_FLAP_M)<=0:
        raise ValueError("Mass, chord and reference quantities must be positive.")
    if not aero_db.mach_bounds[0] <= mach <= aero_db.mach_bounds[1]:
        raise ValueError(f"Mach {mach} outside database bounds {aero_db.mach_bounds}.")
    if aero_db.use_3d_alt and not aero_db.alt_bounds[0] <= alt_m <= aero_db.alt_bounds[1]:
        raise ValueError(f"Altitude {alt_m} outside database bounds {aero_db.alt_bounds}.")
    for pair in ('top','bot'):
        flap_pair_raw(pair,0.,bias,mach,case)
    return lo, hi, bias, sgn, vel, q_inf


def _joint_state(x, mach, alt_m, aero_db, case, q_inf, g_mag, sgn, flight_path_angle_deg):
    state=evaluate_state(float(x[0]),float(x[1]),float(x[2]),mach,alt_m,aero_db,case)
    state.update(acceleration_metrics(state['CL_tot'],state['CD_tot'],q_inf,
                                      g_mag,sgn,flight_path_angle_deg))
    top=state['delta_top_deg']>DELTA_BIAS_OUT_DEG+1e-5
    bot=state['delta_bot_deg']>DELTA_BIAS_OUT_DEG+1e-5
    state['active_pair']='BOTH' if top and bot else ('TOP' if top else ('BOTTOM' if bot else 'NONE'))
    return state


def _joint_seeds(lo, hi, bias, limit, state_at, target):
    """Coarse seeds plus distributed starts; no claim of global optimality."""
    coarse=[]
    for a in np.linspace(lo,hi,7):
        for dt in np.linspace(bias,limit,5):
            for db in np.linspace(bias,limit,5):
                try:
                    state=state_at([a,dt,db])
                except ValueError:
                    continue
                score=(state['Cm_total']/.05)**2+(min(state['requirement_margin_g'],0)/max(target,1))**2
                coarse.append((score,[float(a),float(dt),float(db)]))
    starts=[x for _,x in sorted(coarse,key=lambda item:item[0])[:10]]
    for a in np.linspace(lo,hi,3):
        for dt,db in [(bias,bias),(bias,.5*(bias+limit)),(.5*(bias+limit),bias),
                      (.5*(bias+limit),.5*(bias+limit)),(limit,.5*(bias+limit)),
                      (.5*(bias+limit),limit),(limit,limit)]:
            starts.append([float(a),float(dt),float(db)])
    return list(dict.fromkeys(tuple(x) for x in starts))


def estimate_trim_capability(direction, mach, alt_m, aero_db,
                             max_deflection_deg=None, flight_path_angle_deg=FLIGHT_PATH_ANGLE_DEG,
                             case=NOMINAL_CASE, allocation="both_pairs",
                             max_hinge_torque_Nm=None, initial_guesses=None):
    """Largest normal kinematic turning acceleration found.

    Explicit max_deflection_deg is a CAPABILITY STUDY bound, not a constraint
    imposed on the normal target solver. Hardware limits remain advisory there.
    For one_pair, separately hold each inactive pair at bias and keep the best.
    Nonlinear multistart search; the returned capacity is not a certified maximum.
    """
    lo,hi,bias,sgn,vel,q_inf=_joint_context(1.,direction,mach,alt_m,aero_db,flight_path_angle_deg,case)
    limit=DELTA_MODEL_MAX_DEG if max_deflection_deg is None else float(max_deflection_deg)
    if not bias<limit<=DELTA_MODEL_MAX_DEG:
        raise ValueError('Capability limit must exceed bias and not exceed model search limit.')
    if allocation not in ('both_pairs','one_pair'):
        raise ValueError('allocation must be both_pairs or one_pair.')
    if max_hinge_torque_Nm is not None and (not np.isfinite(max_hinge_torque_Nm) or max_hinge_torque_Nm<=0):
        raise ValueError('Hinge-torque limit must be positive and finite.')
    @lru_cache(maxsize=4096)
    def cached_state(a,dt,db):
        r=_joint_state((a,dt,db),mach,alt_m,aero_db,case,q_inf,1.,sgn,flight_path_angle_deg)
        for pair,delta in (('top',dt),('bot',db)):
            raw=flap_pair_raw(pair,a,delta,mach,case)
            r[f'{pair}_pressure']=raw
            r[f'hinge_torque_{pair}_each_Nm']=abs(q_inf*raw['area_each_m2']*raw['cp_net']*raw['cp_distance_m'])
        return r
    def state_at(x):
        return cached_state(*(float(v) for v in x))
    starts=_joint_seeds(lo,hi,bias,limit,state_at,1.)
    warm=[] if initial_guesses is None else [tuple(np.clip(np.asarray(x,dtype=float),[lo,bias,bias],[hi,limit,limit])) for x in initial_guesses]
    starts=list(dict.fromkeys(warm+starts))
    constraints=[{'type':'eq','fun':lambda x:state_at(x)['Cm_total']/.05},
                 {'type':'ineq','fun':lambda x:state_at(x)['CD_tot']}]
    if max_hinge_torque_Nm is not None:
        constraints.append({'type':'ineq','fun':lambda x:np.array([
            max_hinge_torque_Nm-state_at(x)['hinge_torque_top_each_Nm'],
            max_hinge_torque_Nm-state_at(x)['hinge_torque_bot_each_Nm']])/max_hinge_torque_Nm})
    candidates=[]
    modes=('both',) if allocation=='both_pairs' else ('top','bot')
    for mode in modes:
        bounds=[(lo,hi),(bias,limit),(bias,limit)]
        if mode=='top':bounds[2]=(bias,bias)
        if mode=='bot':bounds[1]=(bias,bias)
        for seed in starts:
            x0=np.array(seed,dtype=float)
            if mode=='top':x0[2]=bias
            if mode=='bot':x0[1]=bias
            try:
                fit=minimize(lambda x:-(state_at(x)['requirement_margin_g']+1),x0,
                             method='SLSQP',bounds=bounds,
                             constraints=constraints,
                             options={'maxiter':180,'ftol':1e-9})
                state=state_at(fit.x).copy()
            except ValueError:
                continue
            if abs(state['Cm_total'])>1e-7 or state['CD_tot'] < -1e-9:
                continue
            if max_hinge_torque_Nm is not None and max(state['hinge_torque_top_each_Nm'],state['hinge_torque_bot_each_Nm'])>max_hinge_torque_Nm+1e-5:
                continue
            capacity=state['requirement_margin_g']+1
            if capacity < -1e-5:continue
            state.update({'capacity_g':max(0.,float(capacity)), 'allocation':allocation,
                          'max_allowed_deflection_deg':limit,'search_method':'multistart_capacity_SLSQP',
                          'solver_nfev':fit.nfev})
            candidates.append(state)
    if not candidates:
        raise TrimSearchError('No verified moment-balanced capability state found; finite search.')
    return max(candidates,key=lambda r:r['capacity_g'])


def _finish_joint_result(result, candidates, g_mag, direction, mach, alt_m,
                         flight_path_angle_deg, case, q_inf, vel, diagnostics):
    result=result.copy()
    sgn=0 if direction.lower()=="either" else (1 if direction.lower()=="nose_up" else -1)
    warnings = ["Quasi-steady estimate; local flow and pressure distribution need calibration."]
    if flight_path_angle_deg is None:
        warnings.append("Flight-path angle was free: reported acceleration uses the selected favorable orientation, not all orientations.")
    if max(result["delta_top_deg"], result["delta_bot_deg"]) > DELTA_ACTUATOR_MAX_DEG:
        warnings.append(f"Required deployment exceeds advisory {DELTA_ACTUATOR_MAX_DEG:g} deg hardware travel.")
    for pair in ("top", "bot"):
        raw = flap_pair_raw(pair, result["alpha_deg"], result[f"delta_{pair}_deg"], mach, case)
        closed = flap_pair_raw(pair, result["alpha_deg"], 0., mach, case)
        result[f"{pair}_pressure"] = raw
        # Gross hinge load on ONE deployed flap, not the incremental vehicle load.
        result[f"hinge_torque_{pair}_each_Nm"] = abs(q_inf*raw["area_each_m2"]*raw["cp_net"]*raw["cp_distance_m"])
        regimes = {raw["regime"], closed["regime"]} if CSV_INCLUDES_CLOSED_OML else {raw["regime"]}
        if "detached_fallback" in regimes:
            warnings.append(f"{pair}: detached-compression modified-Newtonian fallback; pressure fidelity uncertain.")
        if abs(raw['theta_deg']-raw['attached_theta_max_deg']) < 1.0:
            warnings.append(f"{pair}: incidence is within 1 deg of attached-shock limit; model-switch sensitivity.")
        if "vacuum_limit" in regimes:
            warnings.append(f"{pair}: ideal expansion reached vacuum limit.")
        if result[f"delta_{pair}_deg"] >= DELTA_MODEL_MAX_DEG-.1:
            warnings.append(f"{pair}: solution near aerodynamic search boundary.")
    if flight_path_angle_deg is not None and sgn*result["n_normal_g"] < 0:
        warnings.append("Turning direction opposes the requested direction.")
    result.update({"mach": mach, "alt_m": alt_m, "vel_ms": vel,
                   "q_inf_kPa": q_inf/1000., "target_g": g_mag,
                   "flight_path_angle_deg": result.get("effective_flight_path_angle_deg",flight_path_angle_deg),
                   "requested_flight_path_angle_deg":flight_path_angle_deg,
                   "pressure_case": case.name, "warnings": warnings,
                   "search_diagnostics": diagnostics,
                   "saturated": max(result["delta_top_deg"], result["delta_bot_deg"]) > DELTA_ACTUATOR_MAX_DEG,
                   "alpha_extrapolated": False,
                   "L_over_D": abs(result["CL_tot"])/result["CD_tot"] if result["CD_tot"] > 0 else None,
                   "candidate_count": len(candidates),
                   "candidate_solutions": [{k: c[k] for k in ("alpha_deg", "delta_top_deg", "delta_bot_deg", "Cm_total")}
                                           for c in candidates]})
    result["allocation"]="both_pairs"
    result["max_deflection_deg"]=max(result["delta_top_deg"],result["delta_bot_deg"])
    return result

def solve_lifting_body_trim(g_mag, direction, mach, alt_m, aero_db,
                            flight_path_angle_deg=FLIGHT_PATH_ANGLE_DEG,
                            case=NOMINAL_CASE, allocation=CONTROL_ALLOCATION):
    """Find the smallest maximum top/bottom deployment for static target trim.

    Enforce Cm=0 and directional turning acceleration >= target. Minimize
    max(delta_top, delta_bot), then rank ties by their sum. one_pair keeps
    the other pair at its configured bias. Hardware travel stays advisory.
    Search both pressure branches and single-pair boundaries explicitly;
    a finite search gives a best found minimum, not a global certificate.
    """
    if allocation not in ('one_pair','both_pairs'):
        raise ValueError('allocation must be both_pairs or one_pair.')
    lo,hi,bias,sgn,vel,q_inf=_joint_context(g_mag,direction,mach,alt_m,aero_db,flight_path_angle_deg,case)
    limit=DELTA_MODEL_MAX_DEG
    @lru_cache(maxsize=8192)
    def cached_state(a,dt,db):
        return _joint_state((a,dt,db),mach,alt_m,aero_db,case,q_inf,g_mag,sgn,flight_path_angle_deg)
    def state_at(x):
        return cached_state(*(float(v) for v in x)).copy()
    def acceptable(r):
        return (abs(r['Cm_total'])<=1e-7 and r['requirement_margin_g']>=-1e-5
                and r['CD_tot']>=-1e-9)
    candidates=[]
    diagnostics={'objective':'minimize maximum deployment; sum breaks ties',
                 'global_minimum_certified':False}
    def keep(x,method,nfev=0,refined=False):
        x=np.asarray(x,dtype=float)
        if not (lo<=x[0]<=hi and np.all(x[1:]>=bias) and np.all(x[1:]<=limit)):
            return
        r=state_at(x)
        if acceptable(r):
            r.update(search_method=method,solver_nfev=nfev,
                     deployment_minimization_converged=bool(refined))
            candidates.append(r)
    starts=_joint_seeds(lo,hi,bias,limit,state_at,g_mag)
    # Independent single-pair searches also test whether the opposite pair
    # can stay closed. Solve in reduced coordinates, not zero-width bounds.
    modes=('both','top','bot') if allocation=='both_pairs' else ('top','bot')
    for mode in modes:
        def expand(z):
            if mode=='both':return np.asarray(z)
            return np.array([z[0],z[1] if mode=='top' else bias,
                             z[1] if mode=='bot' else bias])
        def residual(z):
            r=state_at(expand(z))
            return [r['Cm_total']/.05,min(r['requirement_margin_g'],0)/max(g_mag,1)]
        seeds=starts if mode=='both' else list(dict.fromkeys((x[0],x[1 if mode=='top' else 2]) for x in starts))
        bounds=([lo,bias,bias],[hi,limit,limit]) if mode=='both' else ([lo,bias],[hi,limit])
        for seed in seeds:
            fit=least_squares(residual,seed,bounds=bounds,x_scale='jac',
                              ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=180)
            keep(expand(fit.x),'static_feasible_'+mode,fit.nfev)
    # The weak-shock / detached fallback switch is discontinuous. A local
    # optimizer can stall there. Parameterize each side of that boundary
    # and solve trim + exact target without treating the pressure jump as a root.
    if allocation=='both_pairs' and case.compression=='shock_expansion':
        for pair in ('top','bot'):
            flow=LOCAL_FLOW[pair]
            me=mach if flow['mach'] is None else float(flow['mach'])
            phi=np.deg2rad(PHI_TOP_DEG if pair=='top' else PHI_BOT_DEG)
            sign=-1 if pair=='top' else 1
            for side in (-1,1):
                theta=attached_shock_limit(me)[1]+side*np.deg2rad(1e-5)
                def expand_boundary(z):
                    a=np.deg2rad(z[0]+flow['alpha_offset_deg']+case.flow_angle_offset_deg)
                    A=np.cos(a);B=sign*np.sin(a)*np.cos(phi)
                    ratio=np.sin(theta)/np.hypot(A,B)
                    if abs(ratio)>1:return None
                    d=np.rad2deg(np.arcsin(ratio)-np.arctan2(B,A))
                    return np.array([z[0],d,z[1]]) if pair=='top' else np.array([z[0],z[1],d])
                def residual_boundary(z):
                    x=expand_boundary(z)
                    if x is None or min(x[1:])<bias or max(x[1:])>limit:
                        return [100.,100.]
                    r=state_at(x)
                    return [r['Cm_total']/.05,r['requirement_margin_g']/max(g_mag,1)]
                for a in np.linspace(lo,hi,5):
                    for d in (bias+.65*(limit-bias),bias+.95*(limit-bias)):
                        fit=least_squares(residual_boundary,[a,d],bounds=([lo,bias],[hi,limit]),
                                          max_nfev=180,ftol=1e-11,xtol=1e-11,gtol=1e-11)
                        x=expand_boundary(fit.x)
                        if x is not None:
                            keep(x,'static_shock_boundary',fit.nfev)
                            if acceptable(state_at(x)):
                                # Refine in boundary coordinates, avoiding the
                                # ill-conditioned three-variable shock-limit gradient.
                                base=np.array([lo,bias]);span=np.array([hi-lo,limit-bias])
                                def bx(y):return expand_boundary(base+span*y)
                                def br(y):return residual_boundary(base+span*y)
                                def delta_bounds(y):
                                    z=bx(y)
                                    return np.r_[z[1:]-bias,limit-z[1:]]/(limit-bias)
                                opt=minimize(lambda y:max(bx(y)[1:])/limit,
                                             (fit.x-base)/span,method='SLSQP',bounds=[(0.,1.)]*2,
                                             constraints=[{'type':'eq','fun':lambda y:br(y)[0]},
                                                          {'type':'ineq','fun':lambda y:br(y)[1]},
                                                          {'type':'ineq','fun':delta_bounds}],
                                             options={'maxiter':150,'ftol':1e-11})
                                keep(bx(opt.x),'static_boundary_minimum_deployment',opt.nfev,opt.success)
    if not candidates:
        # A separate capability search can recover a seed missed by feasibility fitting.
        try:
            cap=estimate_trim_capability(direction,mach,alt_m,aero_db,limit,
                                         flight_path_angle_deg,case,allocation)
            diagnostics['best_capability_found']=cap
            keep([cap['alpha_deg'],cap['delta_top_deg'],cap['delta_bot_deg']],
                 'static_capability_seed',cap['solver_nfev'])
        except TrimSearchError as exc:
            diagnostics['capability_error']=str(exc)
    if not candidates:
        raise TrimSearchError(f'No verified static solution found for {g_mag:g} g with {allocation}. '
                              'Hardware travel was not enforced; finite search is not proof of infeasibility.',diagnostics)
    # Scale all optimization coordinates to [0,1], including the epigraph
    # t=max(top,bottom). Keep each refinement on its seed's pressure branch.
    origin=np.array([lo,bias,bias]);scale=np.array([hi-lo,limit-bias,limit-bias])
    ranked=sorted(candidates,key=lambda r:max(r['delta_top_deg'],r['delta_bot_deg']))
    seeds=[]
    for r in ranked:
        x=np.array([r['alpha_deg'],r['delta_top_deg'],r['delta_bot_deg']])
        if all(np.linalg.norm(x-z)>1. for z in seeds):seeds.append(x)
        if len(seeds)>=12:break
    for x in seeds:
        y=(x-origin)/scale
        opt_bounds=[(0.,1.)]*4
        if allocation=='one_pair':
            opt_bounds[2 if x[2]<=bias+1e-5 else 1]=(0.,0.)
        def physical(z):return origin+scale*z[:3]
        def branch_margins(z):
            values=[]
            for pair,idx in (('top',1),('bot',2)):
                raw=flap_pair_raw(pair,x[0],x[idx],mach,case)
                current=physical(z)
                trial=flap_pair_raw(pair,current[0],current[idx],mach,case)
                side=1 if raw['theta_deg']>raw['attached_theta_max_deg'] else -1
                values.append(side*(trial['theta_deg']-trial['attached_theta_max_deg'])/90.)
            return np.array(values)
        constraints=[{'type':'eq','fun':lambda z:state_at(physical(z))['Cm_total']/.05},
                     {'type':'ineq','fun':lambda z:state_at(physical(z))['requirement_margin_g']/max(g_mag,1)},
                     {'type':'ineq','fun':lambda z:state_at(physical(z))['CD_tot']},
                     {'type':'ineq','fun':lambda z:np.array([z[3]-z[1],z[3]-z[2]])}]
        if case.compression=='shock_expansion':
            constraints.append({'type':'ineq','fun':branch_margins})
        fit=minimize(lambda z:z[3],np.r_[y,max(y[1:])],method='SLSQP',
                     bounds=opt_bounds,constraints=constraints,
                     options={'maxiter':300,'ftol':1e-11})
        # Verify actual physical residuals regardless of optimizer termination.
        keep(physical(fit.x),'static_minimum_deployment',fit.nfev,fit.success)
    # Treat differences below numerical deployment resolution as ties;
    # otherwise a tiny constraint-tolerance undershoot can outrank a
    # successfully minimized state and mislabel the final search status.
    deployment_tolerance_deg=1e-5
    best_max=min(max(r['delta_top_deg'],r['delta_bot_deg']) for r in candidates)
    tied=[r for r in candidates if max(r['delta_top_deg'],r['delta_bot_deg'])<=best_max+deployment_tolerance_deg]
    best_sum=min(r['delta_top_deg']+r['delta_bot_deg'] for r in tied)
    tied=[r for r in tied if r['delta_top_deg']+r['delta_bot_deg']<=best_sum+2*deployment_tolerance_deg]
    result=min(tied,key=lambda r:(not r['deployment_minimization_converged'],abs(r['alpha_deg'])))
    diagnostics['deployment_comparison_tolerance_deg']=deployment_tolerance_deg
    result=_finish_joint_result(result,candidates,g_mag,direction,mach,alt_m,
                                flight_path_angle_deg,case,q_inf,vel,diagnostics)
    result['allocation']=allocation
    result['deployment_objective']='minimize max(top, bottom); sum breaks ties'
    result['minimum_deployment_found_deg']=result['max_deflection_deg']
    result['sum_deployment_deg']=result['delta_top_deg']+result['delta_bot_deg']
    result['deployment_search_status']=('local minimization converged' if result['deployment_minimization_converged']
                                        else 'best feasible candidate; minimization not confirmed')
    return result



def maximize_acceleration(direction, mach, alt_m, aero_db,
                          max_deflection_deg=DELTA_ACTUATOR_MAX_DEG,
                          flight_path_angle_deg=FLIGHT_PATH_ANGLE_DEG,
                          case=NOMINAL_CASE, allocation=CONTROL_ALLOCATION):
    """Maximize normal kinematic turning acceleration, not total acceleration.

    Geometry, mass and flight condition stay fixed. Unlike the target solver,
    this study explicitly bounds deployment by max_deflection_deg (default 45).
    Uniform nominal pressure is used unless an explicit custom case is passed.
    """
    best=estimate_trim_capability(direction,mach,alt_m,aero_db,max_deflection_deg,
                                  flight_path_angle_deg,case,allocation)
    capacity=best['capacity_g']
    if capacity<=0:raise TrimSearchError('No positive joint acceleration capability found.')
    _,_,_,sgn,vel,q_inf=_joint_context(capacity,direction,mach,alt_m,aero_db,flight_path_angle_deg,case)
    x=[best['alpha_deg'],best['delta_top_deg'],best['delta_bot_deg']]
    state=_joint_state(x,mach,alt_m,aero_db,case,q_inf,capacity,sgn,flight_path_angle_deg)
    state.update({'solver_nfev':best['solver_nfev'],'search_method':'maximize_turning_acceleration'})
    result=_finish_joint_result(state,[state],capacity,direction,mach,alt_m,
                                flight_path_angle_deg,case,q_inf,vel,
                                {'objective':'maximize normal kinematic turning acceleration',
                                 'explicit_max_deflection_deg':max_deflection_deg})
    result.update({'capacity_g':capacity,'allocation':allocation})
    return result


def run_sensitivity(g_mag, direction, mach, alt_m, aero_db,
                    flight_path_angle_deg=FLIGHT_PATH_ANGLE_DEG, cases=None, nominal_result=None,
                    allocation=CONTROL_ALLOCATION):
    """Re-solve each case independently. Failures are retained, never dropped.

    These one-at-a-time and mixed cases are illustrative, not an exhaustive
    uncertainty envelope. Cases omit uncertainty in baseline body coefficients.
    """
    rows = []
    for case in PRESSURE_CASES if cases is None else cases:
        try:
            result = (nominal_result if case == NOMINAL_CASE and nominal_result is not None else
                      solve_lifting_body_trim(g_mag, direction, mach, alt_m, aero_db,
                                              flight_path_angle_deg, case, allocation=allocation))
            rows.append({"case": case.name, "status": "solved", "result": result})
        except TrimSearchError as exc:
            rows.append({"case": case.name, "status": "no_solution_found", "reason": str(exc),
                         "diagnostics":exc.diagnostics})
        except ValueError as exc:
            rows.append({"case": case.name, "status": "invalid_configuration", "reason": str(exc)})
    return rows


def _pressure_caveats(state):
    flags=[]
    for pair in ('top','bot'):
        raw=state.get(f'{pair}_pressure')
        if raw is None:continue
        if raw['regime']=='detached_fallback':flags.append(pair+':detached_fallback')
        if abs(raw['theta_deg']-raw['attached_theta_max_deg'])<1.:
            flags.append(pair+':near_shock_switch')
    return ';'.join(flags)


def _refine_envelope_target(cap, target_g, direction, mach, alt_m, aero_db,
                            limit, gamma, case, allocation, torque_limit, objective="min_drag"):
    """Reduce drag (or deployment) from a verified >=target seed, retaining limits.

    Failure to improve retains the verified seed, never changes a met point
    into an unverified optimum. This is a local refinement, not global proof.
    """
    lo,hi,bias,sgn,vel,q=_joint_context(target_g,direction,mach,alt_m,aero_db,gamma,case)
    @lru_cache(maxsize=1024)
    def cached(a,dt,db):
        r=_joint_state((a,dt,db),mach,alt_m,aero_db,case,q,target_g,sgn,gamma)
        for pair,d in (('top',dt),('bot',db)):
            raw=flap_pair_raw(pair,a,d,mach,case)
            r[f'{pair}_pressure']=raw
            r[f'hinge_torque_{pair}_each_Nm']=abs(q*raw['area_each_m2']*raw['cp_net']*raw['cp_distance_m'])
        return r
    def state(x):return cached(*(float(v) for v in x[:3]))
    x0=np.array([cap['alpha_deg'],cap['delta_top_deg'],cap['delta_bot_deg'],
                 max(cap['delta_top_deg'],cap['delta_bot_deg'])])
    bounds=[(lo,hi),(bias,limit),(bias,limit),(bias,limit)]
    if allocation=='one_pair':
        if cap['delta_top_deg']<=bias+1e-5:bounds[1]=(bias,bias)
        else:bounds[2]=(bias,bias)
    cons=[{'type':'eq','fun':lambda x:state(x)['Cm_total']/.05},
          {'type':'ineq','fun':lambda x:state(x)['requirement_margin_g']/max(target_g,1)},
          {'type':'ineq','fun':lambda x:state(x)['CD_tot']},
          {'type':'ineq','fun':lambda x:np.array([x[3]-x[1],x[3]-x[2]])}]
    if torque_limit is not None:
        cons.append({'type':'ineq','fun':lambda x:np.array([
            torque_limit-state(x)['hinge_torque_top_each_Nm'],
            torque_limit-state(x)['hinge_torque_bot_each_Nm']])/torque_limit})
    seed=state(x0).copy()
    try:
        cost=(lambda x:state(x)['CD_tot']/.05) if objective=='min_drag' else (lambda x:x[3])
        fit=minimize(cost,x0,method='SLSQP',bounds=bounds,constraints=cons,
                     options={'ftol':1e-10,'maxiter':220})
        r=state(fit.x).copy()
    except ValueError:
        return seed,False
    valid=(abs(r['Cm_total'])<=1e-7 and r['requirement_margin_g']>=-1e-5 and r['CD_tot']>=-1e-9)
    if torque_limit is not None:
        valid=valid and max(r['hinge_torque_top_each_Nm'],r['hinge_torque_bot_each_Nm'])<=torque_limit+1e-5
    improved=(r['CD_tot']<=seed['CD_tot']+1e-8 if objective=='min_drag' else
              max(r['delta_top_deg'],r['delta_bot_deg'])<=x0[3]+1e-5)
    if valid and improved:return r,True
    return seed,False


def map_turn_envelope(aero_db, mach_grid=None, altitude_grid=None, target_g=15.,
                      max_deflection_deg=DELTA_ACTUATOR_MAX_DEG,
                      flight_path_angle_deg=FLIGHT_PATH_ANGLE_DEG, direction='either',
                      max_q_kPa=None, max_hinge_torque_Nm=None,
                      case=NOMINAL_CASE, allocation=CONTROL_ALLOCATION,
                      refine_target=True, progress=None, flight_path_angle_grid=None, target_objective="min_drag",
                      all_flight_path_angles=False):
    """Quasi-steady feasibility screening across Mach-altitude conditions.

    At EVERY point optimize AoA and both symmetric pairs, not just lift at a
    chosen AoA. Keeps Cm=0, no extrapolation, explicit deployment and optional
    torque/q bounds. If gamma is omitted, a favorable orientation is selected
    per point. This is not a guarantee throughout a turn or a sustained flight.
    'below_target_found' and 'no_verified_trim_found' are finite-search results,
    not proofs of physical infeasibility. Refine the grid near a 15g boundary.
    all_flight_path_angles=True instead screens the worst gravity orientation
    for the same force state, using max(0, abs(L)/(mg)-1) for either direction.
    Target refinement minimizes drag by default; min_deflection is optional.
    """
    if all_flight_path_angles and (flight_path_angle_deg is not None or flight_path_angle_grid is not None):
        raise ValueError('All-angle screening cannot be combined with fixed or sampled flight-path angles.')
    if flight_path_angle_grid is not None:
        if flight_path_angle_deg is not None:
            raise ValueError('Use a fixed flight-path angle or an angle grid, not both.')
        angles=np.unique(np.asarray(flight_path_angle_grid,dtype=float))
        if angles.ndim!=1 or len(angles)==0 or not np.isfinite(angles).all():
            raise ValueError('Flight-path angle grid must be a nonempty finite list.')
        frames=[map_turn_envelope(aero_db,mach_grid,altitude_grid,target_g,max_deflection_deg,
                  float(angle),direction,max_q_kPa,max_hinge_torque_Nm,case,allocation,
                  refine_target,progress,target_objective=target_objective) for angle in angles]
        return pd.concat(frames,ignore_index=True)
    if target_objective not in ('min_drag','min_deflection'):raise ValueError('Unknown target refinement objective.')
    if not np.isfinite(target_g) or target_g<=0:raise ValueError('Target turning g must be positive.')
    if max_q_kPa is not None and (not np.isfinite(max_q_kPa) or max_q_kPa<=0):
        raise ValueError('Maximum dynamic pressure must be positive and finite.')
    if max_hinge_torque_Nm is not None and (not np.isfinite(max_hinge_torque_Nm) or max_hinge_torque_Nm<=0):
        raise ValueError('Maximum hinge torque must be positive and finite.')
    if not DELTA_BIAS_OUT_DEG<max_deflection_deg<=DELTA_MODEL_MAX_DEG:
        raise ValueError('Envelope deployment bound must exceed bias and not exceed model limit.')
    if allocation not in ('both_pairs','one_pair'):raise ValueError('Unknown allocation.')
    if direction not in ('either','nose_up','nose_down'):raise ValueError('Unknown direction.')
    if flight_path_angle_deg is not None and not np.isfinite(flight_path_angle_deg):
        raise ValueError('Flight-path angle must be finite or None.')
    machs=np.unique(np.asarray(mach_grid if mach_grid is not None else
        aero_db.df[aero_db.col_map['mach']].unique(),dtype=float))
    if altitude_grid is None:
        alts=(aero_db.df[aero_db.col_map['alt_m']].unique() if aero_db.use_3d_alt else [ALTITUDE_M])
    else:alts=altitude_grid
    alts=np.unique(np.asarray(alts,dtype=float))
    if machs.ndim!=1 or alts.ndim!=1 or len(machs)==0 or len(alts)==0 or not np.isfinite(machs).all() or not np.isfinite(alts).all():
        raise ValueError('Mach and altitude grids must be nonempty one-dimensional lists.')
    for m in machs:
        validate_freestream_mach(m)
        if not aero_db.mach_bounds[0]<=m<=aero_db.mach_bounds[1]:raise ValueError('Mach grid exceeds database coverage.')
    for h in alts:
        get_atmosphere(h)
        if aero_db.use_3d_alt and not aero_db.alt_bounds[0]<=h<=aero_db.alt_bounds[1]:
            raise ValueError('Altitude grid exceeds database coverage.')
    rows=[]
    for m in machs:
        warm=[]
        for h in alts:
            rho,_,_,sound=get_atmosphere(h);v=m*sound;q=.5*rho*v*v
            row={'Mach':float(m),'altitude_m':float(h),'q_kPa':float(q/1000),
                 'target_turn_g':target_g,'deployment_limit_deg':max_deflection_deg,
                 'hinge_torque_limit_Nm':max_hinge_torque_Nm,'max_q_kPa':max_q_kPa,
                 'allocation':allocation,'direction':direction,'pressure_case':case.name,
                 'target_refinement_objective':target_objective,
                 'flight_path_angle_free':flight_path_angle_deg is None and not all_flight_path_angles,
                 'all_flight_path_angles':all_flight_path_angles,
                 'requested_gamma_deg':flight_path_angle_deg,
                 'mass_kg':MASS_KG,'chord_m':C_FLAP_M,'S_ref_m2':S_REF_M2,
                 'C_ref_m':C_REF_M,'hinge_x_m':X_HINGE_M,
                 'top_area_each_m2':S_FLAP_TOP_EACH_M2,'bottom_area_each_m2':S_FLAP_BOT_EACH_M2,
                 'CG_m':X_CG_M,'phi_top_deg':PHI_TOP_DEG,'phi_bottom_deg':PHI_BOT_DEG,
                 'meets_target':False}
            if max_q_kPa is not None and q/1000>max_q_kPa:
                row['status']='excluded_by_q_limit'
                rows.append(row)
                if progress:progress(row)
                continue
            try:
                cap=estimate_trim_capability(direction,float(m),float(h),aero_db,max_deflection_deg,
                    flight_path_angle_deg,case,allocation,max_hinge_torque_Nm,warm)
            except TrimSearchError as exc:
                row.update(status='no_verified_trim_found',reason=str(exc))
                rows.append(row)
                if progress:progress(row)
                continue
            except ValueError as exc:
                row.update(status='model_error',reason=str(exc))
                rows.append(row)
                if progress:progress(row)
                continue
            warm=[[cap['alpha_deg'],cap['delta_top_deg'],cap['delta_bot_deg']]]
            all_angle_bound=max(0.,abs(cap['n_load_g'])-1.)
            row['same_state_all_angle_turn_lower_bound_g']=all_angle_bound
            row['best_orientation_turn_found_g']=cap['capacity_g'] if flight_path_angle_deg is None else np.nan
            target_gamma=flight_path_angle_deg
            if all_flight_path_angles:
                # Max |L| also maximizes the same-state lower bound over gamma.
                # Evaluate the force state at the least favorable gravity orientation.
                if direction=='nose_up':target_gamma=0.
                elif direction=='nose_down':target_gamma=180.
                else:target_gamma=float(np.rad2deg(np.arccos(np.clip(cap['n_load_g'],-1,1))))
                sgn=0 if direction=='either' else (1 if direction=='nose_up' else -1)
                worst=_joint_state(warm[0],float(m),float(h),aero_db,case,q,1.,sgn,target_gamma)
                cap.update(worst)
                cap['capacity_g']=max(0.,worst['turning_capacity_g'])
            row.update({'status':'target_met' if cap['capacity_g']>=target_g-1e-5 else 'below_target_found',
                        'meets_target':bool(cap['capacity_g']>=target_g-1e-5),
                        'max_turn_found_g':cap['capacity_g'],'margin_g':cap['capacity_g']-target_g,
                        'capacity_alpha_deg':cap['alpha_deg'],
                        'capacity_top_deg':cap['delta_top_deg'],'capacity_bottom_deg':cap['delta_bot_deg'],
                        'capacity_gamma_deg':cap['effective_flight_path_angle_deg'],
                        'capacity_net_accel_g':cap['n_kinematic_g'],'capacity_lift_g':cap['n_load_g'],
                        'capacity_Cm_residual':cap['Cm_total'],
                        'capacity_hinge_top_each_Nm':cap['hinge_torque_top_each_Nm'],
                        'capacity_hinge_bottom_each_Nm':cap['hinge_torque_bot_each_Nm'],
                        'capacity_pressure_caveats':_pressure_caveats(cap)})
            if row['meets_target']:
                target,refined=(_refine_envelope_target(cap,target_g,direction,float(m),float(h),aero_db,
                    max_deflection_deg,target_gamma,case,allocation,max_hinge_torque_Nm,target_objective)
                    if refine_target else (cap,False))
                if all_flight_path_angles and abs(target['n_load_g'])-1 < target_g-1e-5:
                    target,refined=cap,False
                row.update({'target_refined':refined,'target_actual_turn_g':target['turning_capacity_g'],
                            'target_alpha_deg':target['alpha_deg'],
                            'target_top_deg':target['delta_top_deg'],'target_bottom_deg':target['delta_bot_deg'],
                            'target_gamma_deg':target['effective_flight_path_angle_deg'],
                            'target_all_angle_turn_lower_bound_g':max(0.,abs(target['n_load_g'])-1.),
                            'target_Cm_residual':target['Cm_total'],
                            'target_hinge_top_each_Nm':target['hinge_torque_top_each_Nm'],
                            'target_hinge_bottom_each_Nm':target['hinge_torque_bot_each_Nm'],
                            'target_pressure_caveats':_pressure_caveats(target),
                            'target_tangential_g':target['n_tangential_g'],
                            'target_turn_rate_deg_s':float(np.rad2deg(target['n_normal_g']*G_ACCEL/v)),
                            'target_speed_decay_timescale_s':v/(-target['n_tangential_g']*G_ACCEL) if target['n_tangential_g']<0 else np.inf,
                            'target_local_speed_loss_fraction_per_10deg':max(0.,-target['n_tangential_g'])/max(abs(target['n_normal_g']),1e-12)*np.deg2rad(10),
                            'target_drag_to_lift':target['CD_tot']/abs(target['CL_tot']) if abs(target['CL_tot'])>1e-12 else np.inf})
            rows.append(row)
            if progress:progress(row)
    return pd.DataFrame(rows)


def print_turn_envelope(frame):
    print('\nBEST QUASI-STEADY TURN CAPABILITY FOUND [g] — NOT A SUSTAINED-TURN CERTIFICATION')
    if 'max_turn_found_g' in frame:
        for gamma,group in frame.groupby('requested_gamma_deg',dropna=False):
            label=('all orientations, same-state lower bound' if group['all_flight_path_angles'].any() else
                   ('free favorable orientation' if pd.isna(gamma) else f'{gamma:g} deg'))
            print('Flight-path angle: '+label)
            print(group.pivot(index='altitude_m',columns='Mach',values='max_turn_found_g').round(2).to_string())
    print('\nStatus counts: '+str(frame['status'].value_counts().to_dict()))
    met=frame[frame['meets_target']]
    if len(met):
        cols=['Mach','altitude_m','requested_gamma_deg','q_kPa','max_turn_found_g','target_actual_turn_g',
              'target_alpha_deg','target_top_deg','target_bottom_deg','target_gamma_deg',
              'target_hinge_top_each_Nm','target_hinge_bottom_each_Nm','target_speed_decay_timescale_s','target_pressure_caveats']
        print('\nVERIFIED CANDIDATE CONDITIONS MEETING THE TARGET (MODEL AND SUPPLIED LIMITS):')
        print(met[[c for c in cols if c in met]].round(3).to_string(index=False))
    else:print('No verified target-meeting state found on this grid under the supplied limits.')
    if frame['requested_gamma_deg'].nunique()>1:
        intersection=frame.groupby(['Mach','altitude_m'])['meets_target'].all()
        print('\nMACH/ALTITUDE POINTS MEETING THE TARGET AT EVERY SAMPLED FLIGHT-PATH ANGLE:')
        accepted=intersection[intersection]
        print(accepted.to_string() if len(accepted) else 'None verified.')
        print('This tests sampled orientations only; control-rate continuity and duration remain unassessed.')
    print('Refine Mach/altitude sampling around the target boundary; search failures are not proofs of impossibility.')
    if frame['flight_path_angle_free'].any():
        print('Free gamma means a different favorable orientation can be selected at each point.')
    if frame['all_flight_path_angles'].any():
        print('All-angle mode requires the same force state to meet the target for any gravity orientation; Mach/altitude remain fixed at each point.')
    print('Pressure-switch/fallback flags require aerodynamic validation. Duration, rate effects, thermal and structural limits are not assessed.')


def format_net_acceleration_g(result):
    """Table display: actual magnitude if gamma is supplied, otherwise bounds."""
    if result['n_kinematic_g'] is not None:
        return f"{result['n_kinematic_g']:.3f}"
    return f"[{result['n_kinematic_lower_bound_g']:.3f},{result['n_kinematic_upper_bound_g']:.3f}]"


def print_result(res, verbose=False):
    if not verbose:
        print(f"STATIC TRIM — Mach {res['mach']:g}, altitude {res['alt_m']:g} m")
        print(f"  TOP deflection:    {res['delta_top_deg']:.2f} deg")
        print(f"  BOTTOM deflection: {res['delta_bot_deg']:.2f} deg")
        print(f"  AoA: {res['alpha_deg']:+.2f} deg")
        print(f"  Turning acceleration: {res['turning_magnitude_g']:.2f} g "
              f"(requirement/capability level: {res['target_g']:.2f} g)")
        selected=" (freely selected)" if res.get('flight_path_angle_was_selected') else ""
        gamma=res.get('flight_path_angle_deg')
        if gamma is not None:
            print(f"  Flight-path angle: {gamma:+.2f} deg{selected}")
        if 'minimum_deployment_found_deg' in res:
            if res.get('deployment_minimization_converged'):
                print("  Best deployment found; global minimum not certified.")
            else:
                print("  NOTE: Feasible angles found; deployment minimization not confirmed.")
        warnings=res.get('warnings',[])
        if res.get('saturated'):
            print(f"  WARNING: Deployment exceeds advisory {DELTA_ACTUATOR_MAX_DEG:g} deg hardware limit.")
        if any(any(term in w for term in ('fallback','attached-shock limit','vacuum','search boundary')) for w in warnings):
            print("  NOTE: Pressure-model limits affect this result; use --verbose for details.")
        return
    print("="*78)
    print("PITCH-ONLY QUASI-STEADY TRIM ESTIMATE")
    print(f"  Bottom pair moment: 2 * ({CM_DELTA_BOTTOM_EACH_PER_DEG:+.9g} /deg) * deployment_deg")
    print("  Hinge torque and flap-force estimates remain pressure-based.")
    print(f"  Cm contributions (baseline / top / bottom): "
          f"{res['Cm0_cg']:+.7f} / {res['dCm_top']:+.7f} / {res['dCm_bot']:+.7f}")
    print(f"  Allocation: {res.get('allocation','one_pair')}; active: {res['active_pair']}")
    print(f"  Search: {res.get('search_method','unknown')}")
    if 'minimum_deployment_found_deg' in res:
        print(f"  Deployment objective: {res['deployment_objective']}")
        print(f"  Smallest maximum deployment found: {res['minimum_deployment_found_deg']:.4f} deg")
        print(f"  Top + bottom deployment: {res['sum_deployment_deg']:.4f} deg")
        print(f"  Deployment search: {res['deployment_search_status']}; global minimum not certified")
    print(f"  AoA: {res['alpha_deg']:+.4f} deg (no extrapolation)")
    print(f"  Top / bottom outward deployment: {res['delta_top_deg']:.4f} / {res['delta_bot_deg']:.4f} deg")
    print(f"  Turning requirement/capability level, |a_perp|/g: {res['target_g']:.6f} g")
    print(f"  Lift load factor, L/(mg): {res['n_load_g']:+.6f} g")
    print(f"  TURNING acceleration, |a_perp|/g: {res['turning_magnitude_g']:.6f} g")
    print(f"    Signed normal acceleration: {res['n_normal_g']:+.6f} g")
    radius=res['vel_ms']**2/(abs(res['n_normal_g'])*G_ACCEL) if abs(res['n_normal_g'])>1e-12 else np.inf
    print(f"    Instantaneous turn radius at current speed: {radius/1000:.3f} km")
    if res.get('flight_path_angle_was_selected'):
        print("    Flight-path orientation was selected freely for this capability; it is not an all-angle guarantee.")
    if res['n_kinematic_g'] is None:
        print("  Flight-path angle: unspecified")
        print("  Net kinematic acceleration magnitude: not uniquely determined")
        lower=res['n_kinematic_lower_bound_g']
        upper=res['n_kinematic_upper_bound_g']
        print(f"    Possible |a|/g range: [{lower:.6f}, {upper:.6f}] g over all flight-path angles")
        print(f"    Possible |a| range: [{lower*G_ACCEL:.3f}, {upper*G_ACCEL:.3f}] m/s^2")
        print("    Supply --flight-path-angle DEG for the actual magnitude at that angle.")
    else:
        print(f"  Flight-path angle: {res['flight_path_angle_deg']:+.3f} deg")
        print(f"  Net kinematic acceleration magnitude, |a|/g: {res['n_kinematic_g']:.6f} g")
        print(f"  Net kinematic acceleration magnitude, |a|: {res['n_kinematic_magnitude_m_s2']:.3f} m/s^2")
        print(f"    Normal / tangential acceleration: {res['n_normal_g']:+.6f} / {res['n_tangential_g']:+.6f} g")
    print(f"  Cm residual: {res['Cm_total']:+.3e}; requirement margin: {res['requirement_margin_g']:+.3e} g")
    print(f"  Mach {res['mach']:.3f}, altitude {res['alt_m']:.1f} m, q {res['q_inf_kPa']:.3f} kPa")
    print(f"  CL / CD: {res['CL_tot']:+.6f} / {res['CD_tot']:.6f}; verified candidates: {res['candidate_count']}")
    print(f"  Hinge torque per top / bottom flap: {res['hinge_torque_top_each_Nm']:.2f} / {res['hinge_torque_bot_each_Nm']:.2f} Nm")
    for pair in ('top','bot'):
        pressure=res.get(f'{pair}_pressure')
        if pressure is not None:
            print(f"  {pair.upper()} pressure: {pressure['regime']}, incidence={pressure['theta_deg']:.3f} deg, "
                  f"attached limit={pressure['attached_theta_max_deg']:.3f} deg, CP aft of hinge={pressure['cp_distance_m']*1000:.1f} mm")
    for message in res['warnings']:
        print("  NOTE: "+message)


def main():
    global CM_DELTA_BOTTOM_EACH_PER_DEG
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", default=AERO_DATABASE_FILE)
    parser.add_argument("--verbose", action="store_true", help="Show full configuration and aerodynamic diagnostics")
    parser.add_argument("--cm-delta-bottom-per-deg", type=float,
                        default=CM_DELTA_BOTTOM_EACH_PER_DEG,
                        help="Signed vehicle Cm change per degree for ONE bottom flap, about CG; pair contribution is doubled")
    parser.add_argument("--flight-path-angle", type=float, default=FLIGHT_PATH_ANGLE_DEG,
                        help="Velocity angle above local horizontal [deg]; if omitted, choose favorable orientation for any turn")
    parser.add_argument("--g", type=float, default=TARGET_G_MAGNITUDE)
    parser.add_argument("--direction", choices=["either", "nose_up", "nose_down"], default=MANEUVER_DIRECTION)
    parser.add_argument("--mach", type=float, default=MACH_NUMBER)
    parser.add_argument("--altitude", type=float, default=ALTITUDE_M)
    sensitivity_group=parser.add_mutually_exclusive_group()
    sensitivity_group.add_argument("--no-sensitivity", action="store_true")
    sensitivity_group.add_argument("--sensitivity", action="store_true", help="Opt in to pressure sensitivity cases")
    parser.add_argument("--maximize", action="store_true", help="Maximize normal turning acceleration at this flight condition instead of solving --g")
    parser.add_argument("--max-deflection", type=float, default=DELTA_ACTUATOR_MAX_DEG,
                        help="Explicit deployment bound for --maximize; target solver hardware limit stays advisory")
    parser.add_argument("--allocation", choices=["both_pairs","one_pair"], default=CONTROL_ALLOCATION)
    parser.add_argument("--envelope", action="store_true", help="Map optimized turning capability across Mach/altitude")
    parser.add_argument("--mach-grid", type=float, nargs='+', help="Mach samples; default database Mach grid")
    parser.add_argument("--altitude-grid", type=float, nargs='+', help="Altitude samples [m]; default database altitude grid")
    parser.add_argument("--max-q-kpa", type=float, help="Exclude envelope points above this dynamic pressure")
    parser.add_argument("--max-hinge-torque", type=float, help="Per-flap hinge-torque bound [Nm] for --envelope")
    parser.add_argument("--all-flight-path-angles", action="store_true", help="Require the same force state to provide target turning acceleration for any orientation; quasi-steady only")
    parser.add_argument("--flight-path-angle-grid", type=float, nargs='+', help="Flight-path angles [deg] for --envelope; test candidates across turn orientations")
    parser.add_argument("--target-objective", choices=['min_drag','min_deflection'],default='min_drag',
                        help="Refinement objective for target-meeting envelope points")
    parser.add_argument("--envelope-csv", type=Path, help="Save detailed envelope rows to this CSV")
    parser.add_argument("--capability", action="store_true", help="Compare one-pair and four-surface capability at advisory and model limits")
    args = parser.parse_args()
    CM_DELTA_BOTTOM_EACH_PER_DEG=args.cm_delta_bottom_per_deg
    try:
        bottom_pair_cm_increment(0.)
    except (TypeError,ValueError) as exc:
        parser.error(str(exc))
    try:
        validate_freestream_mach(args.mach)
    except ValueError as exc:
        parser.error(str(exc))
    if args.envelope and (args.maximize or args.capability or args.sensitivity):
        parser.error('--envelope is a separate nominal-pressure study; use it without --maximize, --capability or --sensitivity.')
    if args.all_flight_path_angles and not args.envelope:
        parser.error('--all-flight-path-angles applies to --envelope.')
    if args.all_flight_path_angles and (args.flight_path_angle is not None or args.flight_path_angle_grid is not None):
        parser.error('Use --all-flight-path-angles without a fixed or sampled flight-path angle.')
    if args.flight_path_angle_grid is not None and not args.envelope:
        parser.error('--flight-path-angle-grid applies to --envelope.')
    if args.flight_path_angle_grid is not None and args.flight_path_angle is not None:
        parser.error('Use --flight-path-angle-grid or --flight-path-angle, not both.')
    if args.max_hinge_torque is not None and not args.envelope:
        parser.error('--max-hinge-torque currently applies to --envelope.')
    if args.max_q_kpa is not None and not args.envelope:
        parser.error('--max-q-kpa currently applies to --envelope.')
    if args.maximize and args.sensitivity:
        parser.error("--maximize uses the configured nominal pressure case; run --sensitivity separately for a target load.")
    db_path = Path(args.database)
    if not db_path.is_file():
        db_path = Path(__file__).resolve().parent/args.database
    if not db_path.is_file():
        parser.error(f"Aerodynamic database not found: {args.database}. Supply --database PATH.")
    db = FlexibleAeroDatabase.from_file(db_path)
    if args.verbose:
        print(f"CONFIG: chord={C_FLAP_M*1000:.1f} mm, hinge={X_HINGE_M:.4f} m, "
              f"areas(top/bottom each)={S_FLAP_TOP_EACH_M2:.6f}/{S_FLAP_BOT_EACH_M2:.6f} m2")
        task="map turning envelope" if args.envelope else ("maximize turning acceleration" if args.maximize else f"target={args.g:g} g")
        if args.envelope:
            print(f"CONFIG: scanning Mach/altitude grids (database axes unless grid arguments supplied), direction={args.direction}")
        else:
            print(f"CONFIG: Mach={args.mach:g}, altitude={args.altitude:g} m, direction={args.direction}, {task}")
        print(f"CONFIG: mass={MASS_KG:g} kg, body area={S_REF_M2:g} m2, moment length={C_REF_M:g} m, CG={X_CG_M:g} m; fixed baseline coefficients")
        print(f"DATABASE: {db_path.resolve()}; selected columns={db.col_map}")
        print(f"BOTTOM MOMENT: linear per-flap derivative={CM_DELTA_BOTTOM_EACH_PER_DEG:+.9g} /deg; "
              f"pair derivative={2*CM_DELTA_BOTTOM_EACH_PER_DEG:+.9g} /deg; CG reference")
        print("MODEL: bottom moment derivative replaces pressure-derived bottom moment; "
              "top moment, lift/drag and hinge torques retain the pressure model.")
        shape="uniform" if NOMINAL_CASE.k_top==0 and NOMINAL_CASE.k_bot==0 else "linear chordwise gradient"
        print(f"PRESSURE: {shape}; k_top={NOMINAL_CASE.k_top:g}, k_bottom={NOMINAL_CASE.k_bot:g}; Mach domain=[3,8]")
    if args.envelope:
        print(f"ENVELOPE: target={args.g:g} turn g, max deployment={args.max_deflection:g} deg, "
              f"max q={args.max_q_kpa} kPa, per-flap torque limit={args.max_hinge_torque} Nm")
        checkpoints=[]
        if args.envelope_csv is not None:
            args.envelope_csv.parent.mkdir(parents=True,exist_ok=True)
        def report(row):
            checkpoints.append(row.copy())
            if args.envelope_csv is not None:
                pd.DataFrame(checkpoints).to_csv(args.envelope_csv,index=False)
            maximum=row.get('max_turn_found_g')
            detail=f"{maximum:.3f} turn g" if maximum is not None else row.get('reason','')
            print(f"  M={row['Mach']:g}, h={row['altitude_m']:g} m, gamma={row['requested_gamma_deg']}: {row['status']} {detail}",flush=True)
        try:
            frame=map_turn_envelope(db,args.mach_grid,args.altitude_grid,args.g,args.max_deflection,
                args.flight_path_angle,args.direction,args.max_q_kpa,args.max_hinge_torque,
                allocation=args.allocation,progress=report,flight_path_angle_grid=args.flight_path_angle_grid,target_objective=args.target_objective,all_flight_path_angles=args.all_flight_path_angles)
        except ValueError as exc:
            parser.error(str(exc))
        print_turn_envelope(frame)
        if args.envelope_csv is not None:
            frame.to_csv(args.envelope_csv,index=False)
            print(f"Detailed envelope saved to {args.envelope_csv.resolve()}")
        return
    if args.maximize:
        try:
            maximum=maximize_acceleration(args.direction,args.mach,args.altitude,db,
                                           args.max_deflection,args.flight_path_angle,
                                           allocation=args.allocation)
            print_result(maximum, verbose=args.verbose)
            print(f"  Deployment bound: {args.max_deflection:g} deg; global maximum not certified.")
            if args.verbose:
                print(f"  Best turning capability found: {maximum['capacity_g']:.6f} g.")
        except ValueError as exc:
            print(f"MAXIMIZATION: {exc}")
            raise SystemExit(2)
        return
    nominal = None
    try:
        nominal = solve_lifting_body_trim(args.g, args.direction, args.mach, args.altitude, db,
                                         args.flight_path_angle, allocation=args.allocation)
        print_result(nominal, verbose=args.verbose)
    except ValueError as exc:
        print(f"NOMINAL: {exc}")
    if (RUN_SENSITIVITY or args.sensitivity) and not args.no_sensitivity:
        rows = run_sensitivity(args.g, args.direction, args.mach, args.altitude, db,
                               args.flight_path_angle, nominal_result=nominal, allocation=args.allocation)
        print("\nILLUSTRATIVE PRESSURE SENSITIVITY — NOT VALIDATED ACCURACY BOUNDS")
        print(f"{'Case':32} {'AoA':>9} {'Top':>9} {'Bottom':>9} {'Turn |a_perp|/g':>17} {'Lift g':>9} {'Net |a|/g or range':>22}  Advisory hardware")
        solved = []
        for row in rows:
            if row["status"] == "solved":
                r = row["result"]
                solved.append(r)
                print(f"{row['case']:32} {r['alpha_deg']:9.3f} {r['delta_top_deg']:9.3f} "
                      f"{r['delta_bot_deg']:9.3f} {r['turning_magnitude_g']:17.3f} {r['n_load_g']:9.3f} {format_net_acceleration_g(r):>22}  "
                      f"{'EXCEEDED' if r['saturated'] else 'within travel'}")
                for warning in r['warnings']:
                    if any(word in warning for word in ('fallback', 'vacuum', 'boundary', 'opposes')):
                        print('  NOTE: '+warning)
            else:
                print(f"{row['case']:32} {row['status'].upper()}: {row['reason']}")
        print(f"Cases solved: {len(solved)}/{len(rows)}. Finite search; failures do not prove impossibility.")
        if solved:
            for key in ("alpha_deg", "delta_top_deg", "delta_bot_deg"):
                print(f"  {key}: [{min(r[key] for r in solved):.3f}, {max(r[key] for r in solved):.3f}]")
            print("Ranges include solved cases only; inspect failures and model-regime warnings.")
    if args.capability:
        print("\nBEST VERIFIED CAPABILITY FOUND — FINITE SEARCH, NOT CERTIFIED MAXIMUM")
        print(f"{'Allocation':12} {'Limit':>8} {'Turn capability g':>17} {'Lift g':>9} {'Net |a|/g or range':>22} {'AoA':>9} {'Top':>9} {'Bottom':>9}")
        for limit in sorted(set([min(DELTA_ACTUATOR_MAX_DEG,DELTA_MODEL_MAX_DEG),DELTA_MODEL_MAX_DEG])):
            for allocation in ('one_pair','both_pairs'):
                try:
                    cap=estimate_trim_capability(args.direction,args.mach,args.altitude,db,limit,
                                                 args.flight_path_angle,allocation=allocation)
                    print(f"{allocation:12} {limit:8.2f} {cap['capacity_g']:16.3f} {cap['n_load_g']:9.3f} "
                          f"{format_net_acceleration_g(cap):>22} {cap['alpha_deg']:9.3f} "
                          f"{cap['delta_top_deg']:9.3f} {cap['delta_bot_deg']:9.3f}")
                except ValueError as exc:
                    print(f"{allocation:12} {limit:8.2f} NO VERIFIED STATE: {exc}")
    if nominal is None:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
