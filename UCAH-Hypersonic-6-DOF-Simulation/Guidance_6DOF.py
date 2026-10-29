import numpy as np
from Atmosphere_6DOF import getAtm
from Vehicle_6DOF import TRAJECTORY
class PitchAutopilot:
    def __init__(self, dt):
        self.dt = dt
        
        # Outer Loop: Altitude error -> Alpha command
        self.K_z = 0.002   
        
        # Inner Loop: Alpha error -> Virtual Pitch Moment
        self.Kp_alpha = 3.5 
        self.Kd_q = 1.5 
        
        self.q_design = 50000.0 
        self.max_deflection = np.deg2rad(20.0)
        self.max_rate = np.deg2rad(200.0)
        
        # State tracking for exactly 2 surfaces: [Deck, Keel]
        self.prev_fins = np.zeros(2) 

    def update(self, state, target_h):
        # 1. Extract pure pitch kinematics from the state array
        x = state[0]
        z = state[2]
        u = state[3]
        w = state[5]
        q = state[11] # Pitch rate
        
        # 2. Compute longitudinal aerodynamic environment
        h = -z
        V = np.hypot(u, w) 
        alpha = np.arctan2(w, u)
        
        # Get density for dynamic pressure gain scaling
        _, _, rho, _ = getAtm(h)
        qbar = 0.5 * rho * V**2
        
        # 3. Outer Loop: Trajectory Controller
        err_h = target_h - h
        alpha_cmd = self.K_z * err_h
        
        # 4. Inner Loop: Pitch Attitude Controller
        q_scale = self.q_design / max(qbar, 1000.0) 
        pitch_cmd = (self.Kp_alpha * (alpha_cmd - alpha) - self.Kd_q * q) * q_scale
        
        # 5. Pure Pitch Mixer (Positive flush deflections only)
        # Pitch Up: Deploy Deck. Pitch Down: Deploy Keel.
        cmd_Deck = max(0.0, pitch_cmd)
        cmd_Keel = max(0.0, -pitch_cmd)
        
        raw_fins = np.array([cmd_Deck, cmd_Keel])
        
        # 6. Actuator Limits
        actual_fins = np.zeros(2)
        for i in range(2):
            # Apply mechanical rate limits
            max_delta = self.max_rate * self.dt
            desired_delta = raw_fins[i] - self.prev_fins[i]
            clipped_delta = np.clip(desired_delta, -max_delta, max_delta)
            
            # Apply absolute deflection limits (Hard floor at 0.0)
            pos = self.prev_fins[i] + clipped_delta
            actual_fins[i] = np.clip(pos, 0.0, self.max_deflection)
            
        self.prev_fins = actual_fins
        return actual_fins

# Initialize globally with your simulation timestep
autopilot = PitchAutopilot(dt=0.01)

def getControls(t, state):
    # Extract downrange distance (x) from the state array[cite: 1]
    downrange_x = state[0]
    
    # Placeholder Trajectory: Replace these with your actual glide path profile
    traj_x = TRAJECTORY[:, 0]
    traj_h = TRAJECTORY[:, 1]
    traj_t = TRAJECTORY[:, 2]
    
    # Interpolate the target altitude for the current downrange position
    target_h = np.interp(downrange_x, traj_x, traj_h)
    
    # Run the pitch autopilot
    deltas = autopilot.update(state, target_h)
    
    # Returns [deck_delta, keel_delta] in radians
    return deltas



if __name__ == "__main__":
    import numpy as np
    
    # Mock getAtm if it isn't imported in your test environment
    # Returns: T, P, rho, a
    def mock_getAtm(h):
        return 216.65, 5474.8, 0.088, 295.0 # Roughly 20km altitude conditions
    
    # Temporarily override the global getAtm for the test
    getAtm

    print("--- Starting Pitch Autopilot Unit Test ---\n")
    
    # 1. Initialize Controller
    dt = 0.01
    test_autopilot = PitchAutopilot(dt=dt)
    
    # 2. Setup Dummy State (14 elements)
    # Missile is at 20,000m altitude (z = -20000)
    # Moving forward at Mach 5 (u = 1500 m/s)
    # Flying perfectly level (w = 0, q = 0)
    dummy_state = np.zeros(14)
    dummy_state[0] = 5000.0   # x: 5km downrange
    dummy_state[2] = -20000.0 # z: 20km altitude
    dummy_state[3] = 1500.0   # u: 1500 m/s forward
    dummy_state[5] = 0.0      # w: 0 m/s vertical
    dummy_state[11] = 0.0     # q: 0 rad/s pitch rate
    
    # 3. Setup Target
    target_altitude = 25000.0 # We want to climb 5,000m
    
    print(f"Current Altitude : {-dummy_state[2]} m")
    print(f"Target Altitude  : {target_altitude} m")
    print(f"Altitude Error   : {target_altitude - (-dummy_state[2])} m -> Expecting Pitch UP (Deck Deployment)\n")
    
    print("Step | Deck Cmd (deg) | Keel Cmd (deg) | Notes")
    print("-" * 55)
    
    # 4. Run loop for 10 timesteps (0.1 seconds of flight)
    for step in range(1, 11):
        # Run autopilot
        fins_rad = test_autopilot.update(dummy_state, target_altitude)
        
        # Convert output to degrees for easy reading
        deck_deg = np.degrees(fins_rad[0])
        keel_deg = np.degrees(fins_rad[1])
        
        # Determine what the limits are doing
        notes = ""
        if step == 1:
            notes = "First step: limited by max_rate (200 deg/s * 0.01s = 2.0 deg)"
        elif step == 10:
            notes = "Hitting max deflection limit (20.0 deg)" if deck_deg >= 20.0 else "Ramping up..."
            
        print(f"{step:4d} | {deck_deg:14.2f} | {keel_deg:14.2f} | {notes}")
        
    print("\nTest Passed If:")
    print("1. Keel stays exactly at 0.00 (Flaps only push outward).")
    print("2. Deck ramps up by exactly 2.0 degrees per step (Rate limits working).")
    print("3. Deck caps at 20.0 degrees (Deflection limit working).")