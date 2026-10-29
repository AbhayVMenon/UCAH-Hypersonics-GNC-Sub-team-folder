# 1. make document of inputs and outputs

# 2. make adaptations to have seperate deck/keel control surfaces

# 3. provode variable diagram for control surface inputs and outputs

### place ohlder for now.
from AnimateTraj import animate_trajectory
from Sim_6DOF import run6DOF
from Atmosphere_6DOF import getAtm
import numpy as np
from display_6DOF import plotAltitude, plotAttitude, plotDiagnostics, plotRun, plotGround, plot3D, plotAOA
from Vehicle_6DOF import vehicleParams
launch_angle = np.deg2rad(25)           # launch flight path angle (rad)
h0 = 1e3                                # launch altitude, m
x0 = 0.0                                # downrange start, m
y0 = 0.0                                # crossrange start, m
T, P, rho, a = getAtm(h0)               # local conditions at launch
v0 = 8 * a                              # Mach 8 launch speed, m/s
simParams = {
    "dt": 0.025,         # integration timestep, s
    "t_max": 600.0,     # 10 min cutoff, s
}

initial_state = [ x0, y0, -h0,                                      # position, flat-Earth NED (m)
          v0, 0, 0,                                                 # velocity in BODY axes (m/s)
          np.cos(launch_angle/2), 0, np.sin(launch_angle/2), 0,     # attitude quaternion, scalar-first
          0, 0, 0,                                                  # body angular rates p, q, r (rad/s)
          0 ]                                                       # actual fin deflection (rad)
history, state, t = run6DOF(initial_state, vehicleParams=vehicleParams, simParams=simParams )
plotAltitude(history)
plotRun(history)
# plotGround(history)
plot3D(history)
plotDiagnostics(history)
plotAttitude(history)
# animate_trajectory(history)
plotAOA(history)