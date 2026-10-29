 # animate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Frame_Math import dcm_from_quat
from Vehicle_6DOF import vehicleParams

def animate_trajectory(history, skip_frames=10):
    """
    Animates the 6DOF trajectory using matplotlib.
    skip_frames: speeds up animation by skipping integration steps.
    """
    # Extract states (NED positions and scalar-first quaternions)
    x_ned = history[::skip_frames, 1]
    y_ned = history[::skip_frames, 2]
    z_ned = history[::skip_frames, 3]
    quats = history[::skip_frames, 6:10]
    
    # 1. Build a basic 3D body in body axes (x out nose, y right, z down)
    # Using your dimensions from Vehicle_6DOF
    L = vehicleParams["length"]
    R = vehicleParams["d_ref"] / 2
    x_cg = vehicleParams["x_cg"]
    
    # Create a simple cylindrical body mesh relative to CG
    X_body = np.linspace(-L + x_cg, x_cg, 10)  # nose is at +x_cg, tail at -L+x_cg
    theta = np.linspace(0, 2*np.pi, 8)
    X, T = np.meshgrid(X_body, theta)
    Y = R * np.cos(T)
    Z = R * np.sin(T)
    
    # Flatten for matrix multiplication
    body_points = np.vstack([X.flatten(), Y.flatten(), Z.flatten()])

    # 2. Set up the plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('North (m)')
    ax.set_ylabel('East (m)')
    ax.set_zlabel('Altitude (m)')
    
    # Set fixed axes limits based on trajectory bounds
    ax.set_xlim([min(x_ned), max(x_ned)])
    ax.set_ylim([min(y_ned)-100, max(y_ned)+100])
    ax.set_zlim([0, max(-z_ned)])
    ax.set_box_aspect([1, 1, 0.5]) # Squashed Z for better viewing
    
    # Plot elements: The full trajectory line, and the moving 3D vehicle wireframe
    ax.plot(x_ned, y_ned, -z_ned, 'k--', alpha=0.5, label="Trajectory")
    vehicle_mesh, = ax.plot([], [], [], color='b', marker='.', linestyle='none', markersize=2)

    # 3. Animation update function
    def update(frame):
        # Current position (invert Z for standard altitude rendering)
        pos_ned = np.array([x_ned[frame], y_ned[frame], z_ned[frame]])
        
        # Get DCM (Body -> NED)
        C_bn = dcm_from_quat(quats[frame])
        
        # Rotate body points to NED frame, then add translation
        earth_points = C_bn @ body_points
        
        # Shift to current Earth position (invert Z for plotting altitude upward)
        plot_x = earth_points[0, :] + pos_ned[0]
        plot_y = earth_points[1, :] + pos_ned[1]
        plot_z = -(earth_points[2, :] + pos_ned[2]) 
        
        vehicle_mesh.set_data(plot_x, plot_y)
        vehicle_mesh.set_3d_properties(plot_z)
        
        return vehicle_mesh,

    ani = animation.FuncAnimation(fig, update, frames=len(x_ned), interval=20, blit=False)
    plt.show()