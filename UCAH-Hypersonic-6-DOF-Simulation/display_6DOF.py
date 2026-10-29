import numpy as np
import matplotlib.pyplot as plt

from Vehicle_6DOF import MASS, G0, S_REF
from Aerodynamics_6DOF import getAero

MACH_REQ = 3.0
CEILING = 30.0                          # km
ALPHA_MAX = np.deg2rad(25.0)            # trimmable limit -- unverified, needs C_M

# q floor for a 15 g pull: 15*m*g = q*S*CL_max
# computed, not hardcoded, so it tracks changes in mass, area, and aero model
CL_MAX, _ = getAero(8.0, ALPHA_MAX)
Q_FLOOR = 15.0 * MASS * G0 / (S_REF * CL_MAX)


def plotRun(history):
    (t, x, y, z, u, v, w, q0, q1, q2, q3, roll_rate, pitch_rate, yaw_rate, delta,
     h, V, mach, qbar, alpha, beta, nz) = history.T
    
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    ax[0,0].plot(np.hypot(x, y)/1000, h/1000)
    ax[0,0].axhline(CEILING, ls='--', c='r', lw=0.8)
    ax[0,0].set_xlabel('ground range, km'); ax[0,0].set_ylabel('altitude, km')
    ax[0,0].set_title('Trajectory')

    ax[0,1].plot(t, mach)
    ax[0,1].axhline(MACH_REQ, ls='--', c='r', lw=0.8)
    ax[0,1].set_xlabel('time, s'); ax[0,1].set_ylabel('Mach')
    ax[0,1].set_title(f'Impact M = {mach[-1]:.2f}')

    ax[1,0].plot(t, qbar/1000)
    #ax[1,0].axhline(Q_FLOOR/1000, ls='--', c='r', lw=0.8)
    ax[1,0].set_xlabel('time, s'); ax[1,0].set_ylabel('q, kPa')
    ax[1,0].set_title(f'Dynamic Pressure)')

    ax[1,1].plot(t, nz)
    ax[1,1].axhline(15, ls='--', c='r', lw=0.8)
    ax[1,1].set_xlabel('time, s'); ax[1,1].set_ylabel('load factor, g')
    ax[1,1].set_title('Normal Acceleration')

    for a_ in ax.flat:
        a_.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


def plotGround(history):
    """Crossrange vs downrange. Flat line until something drives y."""
    (t, x, y, z, u, v, w, q0, q1, q2, q3, roll_rate, pitch_rate, yaw_rate, delta,
     h, V, mach, qbar, alpha, beta, nz) = history.T

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(x/1000, y/1000)
    ax.set_xlabel('downrange x, km'); ax.set_ylabel('crossrange y, km')
    ax.set_title('ground track')
    ax.axis('equal'); ax.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


def plot3D(history):
    """3D trajectory. Sits in the y = y0 plane until something drives y."""
    (t, x, y, z, u, v, w, q0, q1, q2, q3, roll_rate, pitch_rate, yaw_rate, delta,
     h, V, mach, qbar, alpha, beta, nz) = history.T

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection='3d')

    ax.plot(x/1000, y/1000, h/1000, lw=1.5)
    ax.plot(x/1000, y/1000, 0, ls='--', c='gray', lw=0.8)       # ground shadow

    ax.scatter(x[0]/1000, y[0]/1000, h[0]/1000, c='g', s=30, label='launch')
    ax.scatter(x[-1]/1000, y[-1]/1000, 0, c='r', s=30, label='impact')

    ax.set_xlabel('downrange x, km')
    ax.set_ylabel('crossrange y, km')
    ax.set_zlabel('altitude, km')
    ax.set_zlim(0, None)
    ax.legend()
    fig.tight_layout()
    plt.show()

def plotAOA(history):
    """Angle of attack vs time."""
    (t, x, y, z, u, v, w, q0, q1, q2, q3, roll_rate, pitch_rate, yaw_rate, delta,
     h, V, mach, qbar, alpha, beta, nz) = history.T

    a_deg = np.degrees(alpha)                                   # rad -> deg for reading

    # startup transient is ~1000x the settled alpha, so zoom in after it dies out:
    # settled = after the last point where |alpha| exceeds 1.5x its max over the back half of the run
    settled_max = np.abs(a_deg[len(a_deg)//2:]).max()
    big = np.nonzero(np.abs(a_deg) > 1.5 * settled_max)[0]
    i0 = big[-1] + 1 if big.size and big[-1] + 1 < len(t) else 0

    # no ALPHA_MAX lines here: +/-25 deg would squash the tiny alpha variations flat
    fig, ax = plt.subplots(2, 1, figsize=(8, 6))

    ax[0].plot(t, a_deg)
    ax[0].axvspan(t[i0], t[-1], color='gray', alpha=0.1)       # region shown in the zoom below
    ax[0].set_title('Angle of Attack, full run')

    ax[1].plot(t[i0:], a_deg[i0:])
    ax[1].set_title(f'Angle of Attack, t > {t[i0]:.1f} s (settled)')
    ax[1].set_xlabel('time, s')

    for a_ in ax:
        a_.margins(y=0.1)                                       # autoscale tight to the data
        a_.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        a_.set_ylabel('alpha, deg')
        a_.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


def plotAltitude(history):
    t = history[:, 0]
    alt = history[:, 15] / 1000.0  # Extract 'h' and convert to km
    
    plt.figure()
    plt.plot(t, alt, linewidth=2)
    plt.title("Altitude vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (km)")
    plt.grid(True, alpha=0.3)
    plt.show()


import matplotlib.pyplot as plt
import numpy as np

def plotDiagnostics(history):
    """
    Analyzes the flight history to verify aerodynamic responses to fin deflections.
    Assumes history columns:
    0: time
    10, 11, 12: roll rate (p), pitch rate (q), yaw rate (r)
    18, 19: alpha (rad), beta (rad) (Assuming these follow h, V, mach, qbar in outputs)
    """
    t = history[:, 0]
    
    # State vector is columns 1-14. Rates are 11, 12, 13.
    p = np.degrees(history[:, 11])  # Roll rate
    q = np.degrees(history[:, 12])  # Pitch rate
    r = np.degrees(history[:, 13])  # Yaw rate
    
    # Outputs dict starts at col 15. 
    # (15: h, 16: V, 17: mach, 18: qbar, 19: alpha, 20: beta, 21: nz)
    alpha = np.degrees(history[:, 19]) 
    beta = np.degrees(history[:, 20])

    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.suptitle("6DOF Flight Dynamics Validation", fontsize=14, fontweight='bold')

    # Panel 1: Pitch Analysis
    axs[0].plot(t, q, label="Pitch Rate (q) [deg/s]", color='b')
    axs[0].plot(t, alpha, label="Angle of Attack (alpha) [deg]", color='c', linestyle='--')
    axs[0].set_ylabel("Pitch Axis")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(loc="upper left")
    axs[0].axvspan(5, 15, color='gray', alpha=0.1)   # Pitch Up region
    axs[0].axvspan(15, 25, color='gray', alpha=0.2)  # Pitch Down region
    
    # Panel 2: Yaw Analysis
    axs[1].plot(t, r, label="Yaw Rate (r) [deg/s]", color='r')
    axs[1].plot(t, beta, label="Sideslip (beta) [deg]", color='m', linestyle='--')
    axs[1].set_ylabel("Yaw Axis")
    axs[1].grid(True, alpha=0.3)
    axs[1].legend(loc="upper left")
    axs[1].axvspan(25, 35, color='gray', alpha=0.1)  # Yaw Right region
    axs[1].axvspan(35, 45, color='gray', alpha=0.2)  # Yaw Left region

    # Panel 3: Roll Analysis
    axs[2].plot(t, p, label="Roll Rate (p) [deg/s]", color='g')
    axs[2].set_ylabel("Roll Axis")
    axs[2].set_xlabel("Time (s)")
    axs[2].grid(True, alpha=0.3)
    axs[2].legend(loc="upper left")
    axs[2].axvspan(45, t[-1], color='gray', alpha=0.1) # Roll region

    plt.tight_layout()
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def plotAttitude(history):
    """
    Calculates and plots the Euler angles (Pitch, Roll, Yaw) over time 
    from the quaternion state history.
    """
    # Unpack time and quaternion columns from the history array[cite: 5]
    t = history[:, 0]
    q0 = history[:, 7]
    q1 = history[:, 8]
    q2 = history[:, 9]
    q3 = history[:, 10]

    # Convert quaternions to Euler angles (Z-Y-X aerospace sequence, in radians)
    roll = np.arctan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1**2 + q2**2))
    
    # Calculate pitch, clamping the input to arcsin to prevent floating point domain errors
    sinp = 2 * (q0 * q2 - q3 * q1)
    sinp = np.clip(sinp, -1.0, 1.0) 
    pitch = np.arcsin(sinp)
    
    yaw = np.arctan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2**2 + q3**2))

    # Convert radians to degrees for readability
    roll_deg = np.degrees(roll)
    pitch_deg = np.degrees(pitch)
    yaw_deg = np.degrees(yaw)

    # Plot the three axes in stacked subplots
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Vehicle Attitude (Euler Angles)", fontsize=14, fontweight='bold')

    # Pitch Panel
    axs[0].plot(t, pitch_deg, label="Pitch (Theta)", color='b')
    axs[0].set_ylabel("Pitch [deg]")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(loc="upper right")

    # Roll Panel
    axs[1].plot(t, roll_deg, label="Roll (Phi)", color='g')
    axs[1].set_ylabel("Roll [deg]")
    axs[1].grid(True, alpha=0.3)
    axs[1].legend(loc="upper right")

    # Yaw Panel
    axs[2].plot(t, yaw_deg, label="Yaw (Psi)", color='r')
    axs[2].set_ylabel("Yaw [deg]")
    axs[2].set_xlabel("Time (s)")
    axs[2].grid(True, alpha=0.3)
    axs[2].legend(loc="upper right")

    plt.tight_layout()
    plt.show()