import numpy as np
from Vehicle_6DOF import vehicleParams
SHOW_PLOT = True  # True -> open the interactive 3D fin/body debug plot

# ---------- edit this part ----------

vehicleParams["c_g"] = 0.0   # initial fin deflection (rad)

theta = vehicleParams["fin_theta"]              # deg
phi = np.radians(90-theta)        # the angle that stays fixed
# fin 1 hinge relative to the CG (origin): aft (-N), +E, up (-D)
L1 = -vehicleParams["fin_from_cg_n"]
L2 = vehicleParams["fin_from_cg_e"]
L3 = -vehicleParams["fin_from_cg_d"]
h = np.sqrt(vehicleParams["fin_area"])          # square flap: side = sqrt(area)
finMid = vehicleParams["fin_center"]            # hinge -> flap midpoint (= h/2)
nose_len = 1 #vehicleParams["x_cg"] - L1           # nose tip to fin bases (+N)
tail_len = .5 #vehicleParams["length"] - nose_len   # fin bases to aft end (-N)
def position(delta, phi):
    # Replace with your own vector components.
    # Example: arm rotating by delta about z, second link tilted up by phi
    N = (L1 + -1*finMid * np.cos(delta))
    E = (L2 + finMid * np.cos(phi) * np.sin(delta))
    D = L3 + -1 * finMid * np.sin(phi) * np.sin(delta)
    return np.array([N, E, D])
# ------------------------------------

hinge = np.array([L1, L2, L3])

def flap_corners(delta, phi):
    # chord direction: hinge -> flap midpoint
    r = (position(delta, phi) - hinge) / finMid
    # hinge axis: normal to the plane the midpoint sweeps through
    a = np.cross(position(0, phi) - hinge, position(np.pi/2, phi) - hinge)
    a /= np.linalg.norm(a)
    # square flap: chord 0..h along r, span +-h/2 along the hinge axis
    return np.array([hinge + s*h*r + t*(h/2)*a
                     for s, t in [(0, -1), (1, -1), (1, 1), (0, 1)]])

# Mirror sign on (N, E, D) for each fin, numbered like Cartesian quadrants
# when looking forward from behind (E to the right, up = -D):
# fin 1 = +E/up, fin 2 = -E/up, fin 3 = -E/down, fin 4 = +E/down
mirrors = [np.array([1,  1,  1]),
           np.array([1, -1,  1]),
           np.array([1, -1, -1]),
           np.array([1,  1, -1])]

def fin_r(delta, fin, phi=phi):
    # r vector (N, E, D) to the flap midpoint of fin 1-4 (quadrant numbering above)
    # delta in radians
    if fin not in (1, 2, 3, 4):
        raise ValueError("fin must be 1, 2, 3 or 4")
    return position(delta, phi) * mirrors[fin - 1]


def plot_fins():
    # Interactive debug view: body of revolution, 4 fins, and a slider per flap.
    # Only imported/run when SHOW_PLOT is True.
    import matplotlib
    matplotlib.use("TkAgg")      # standalone window that works outside Jupyter/Interactive
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    colors = ["g", "b", "m", "orange"]
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    plt.subplots_adjust(bottom=0.3)

    # Full path of each fin midpoint, for reference
    th_all = np.linspace(0, 0.5*np.pi, 300)
    path = np.array([position(th, phi) for th in th_all]).T

    vecs, dots, flaps = [], [], []
    for m, c in zip(mirrors, colors):
        ax.plot(*(path * m[:, None]), "r-", alpha=0.3)

        # The position vector and its tip
        p = position(0, phi) * m
        vec, = ax.plot([0, p[0]], [0, p[1]], [0, p[2]], "-", color=c, lw=2)
        dot, = ax.plot([p[0]], [p[1]], [p[2]], "ro")

        # The flap surface
        flap = Poly3DCollection([flap_corners(0, phi) * m], facecolor=c, edgecolor="k", alpha=0.4)
        ax.add_collection3d(flap)

        vecs.append(vec); dots.append(dot); flaps.append(flap)

    # Body of revolution about the N axis, radius chosen so it passes through all 4 fin bases.
    # Tangent ogive nose from the tip down to the bases, then a straight cylinder toward -N.
    # Depends only on L1, L2, L3 (and nose_len / tail_len), not on delta.
    R = np.hypot(L2, L3)
    rho = (R**2 + nose_len**2) / (2*R)                     # ogive arc radius
    x_nose = np.linspace(0, nose_len, 60)                  # distance aft of the tip
    r_nose = np.sqrt(rho**2 - (nose_len - x_nose)**2) + R - rho
    N_body = np.concatenate([L1 + nose_len - x_nose, [L1 - tail_len]])
    r_body = np.concatenate([r_nose, [R]])
    ang = np.linspace(0, 2*np.pi, 60)                      # sweep angle about N (don't reuse phi)
    NN, PP = np.meshgrid(N_body, ang)
    RR = np.meshgrid(r_body, ang)[0]
    ax.plot_surface(NN, RR*np.cos(PP), RR*np.sin(PP), color="gray", alpha=0.25, linewidth=0)
    # outline: profile lines along the body, and rings at the fin bases and the aft end
    for ph in np.linspace(0, 2*np.pi, 8, endpoint=False):
        ax.plot(N_body, r_body*np.cos(ph), r_body*np.sin(ph), "k-", lw=0.5, alpha=0.5)
    for Nr in (L1, L1 - tail_len):
        ax.plot(np.full_like(ang, Nr), R*np.cos(ang), R*np.sin(ang), "k-", lw=0.8)

    # Equal scale on all axes, fitted to the body and fins
    pts = np.column_stack([[N_body.max(), R, R], [N_body.min(), -R, -R],
                           *[(flap_corners(d, phi) * m).T for m in mirrors for d in th_all[::30]]])
    lo, hi = pts.min(axis=1) - 0.1, pts.max(axis=1) + 0.1
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(hi - lo)
    ax.set_xlabel("N"); ax.set_ylabel("E"); ax.set_zlabel("D")

    sliders = []
    for i in range(4):
        slider_ax = fig.add_axes([0.2, 0.17 - 0.05*i, 0.6, 0.03])
        sliders.append(Slider(slider_ax, f"δ{i+1} (deg)", 0, 90, valinit=0))

    def make_update(i):
        def update(val):
            delta = np.radians(val)
            m = mirrors[i]
            p = position(delta, phi) * m
            vecs[i].set_data_3d([0, p[0]], [0, p[1]], [0, p[2]])
            dots[i].set_data_3d([p[0]], [p[1]], [p[2]])
            flaps[i].set_verts([flap_corners(delta, phi) * m])
            fig.canvas.draw_idle()
        return update

    for i, slider in enumerate(sliders):
        slider.on_changed(make_update(i))
    fig._sliders = sliders       # keep a reference so the sliders aren't garbage-collected

    plt.show(block=True)


def fin_loads(q, alpha, beta, deltas, mach, gamma=1.4, r_ref=np.zeros(3), rate=None, V=None):
    # q      : dynamic pressure [Pa]
    # alpha  : angle of attack [rad]
    # beta   : sideslip [rad]
    # deltas : 4 flap deflections [rad], fins 1-4 (quadrant numbering above)
    # mach   : None -> Newtonian (Cp_max = 2), otherwise modified Newtonian
    # r_ref  : moment reference point (N, E, D), relative to the CG (fin geometry origin)
    # rate   : body rates (p, q, r) [rad/s]; with V, adds omega x r at each flap (fin damping)
    # V      : airspeed [m/s], needed with rate
    # returns F [N] and M [N*m] as (N, E, D) vectors
    #plot_fins()
    if mach is None:
        cp_max = 2.0
    else:
        g, M2 = gamma, mach**2
        p02_p1 = ((g+1)**2*M2 / (4*g*M2 - 2*(g-1)))**(g/(g-1)) * (1 - g + 2*g*M2) / (g+1)
        cp_max = 2/(g*M2) * (p02_p1 - 1)

    # direction the vehicle moves through the air, in body axes
    v_hat = np.array([np.cos(alpha)*np.cos(beta),
                      np.sin(beta),
                      np.sin(alpha)*np.cos(beta)])
    

    F = np.zeros(3)
    M = np.zeros(3)
    for fin, delta in enumerate(deltas, start=1):
        m = mirrors[fin - 1]
        c = flap_corners(delta, phi) * m
        r_fin = fin_r(delta, fin, phi)
        n = np.cross(c[1] - c[0], c[3] - c[0])   # panel normal, length = area
        area = np.linalg.norm(n)
        n /= area
        # mirroring flips the corner order on fins 2 and 4, so force n onto the exposed
        # (outer) face: it points forward and away from the body for 0 <= delta <= 90 deg
        hinge_m = hinge * m
        if n @ np.array([1.0, hinge_m[1], hinge_m[2]]) < 0:
            n = -n
        # local flow direction at the flap: freestream plus the body rotation (omega x r)
        v_loc = v_hat if rate is None or not V else v_hat + np.cross(rate, r_fin) / V
        s = v_loc @ n                            # sin of flow incidence angle on the outer face
        # Newtonian: Cp = Cp_max*sin^2 on the windward face, 0 when shadowed (s < 0)
        dF = -q * cp_max * area * max(s, 0.0)**2 * n
        F += dF
        M += np.cross(r_fin - r_ref, dF)  # CP at flap midpoint
    return F, M


# q = 50000.0
# alpha = np.radians(10)
# beta = np.radians(10)
# delta = np.radians(np.array([10, 10, 0, 10]))
# mach = 5.0
# F, M = fin_loads(q, alpha, beta, delta, mach=mach)
# print("F = ", F)
# print("M = ", M)

# if SHOW_PLOT:
#     plot_fins()