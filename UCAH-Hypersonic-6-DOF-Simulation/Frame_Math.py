## 
"""
Frame_Math.py -- attitude math for the 6DOF.
Conventions (every function in this file assumes these):
  - Quaternion q = [q0, q1, q2, q3], SCALAR FIRST (q0 = cos(angle/2))
  - q rotates BODY -> NED  (C_bn @ v_body = v_ned)
  - NED: x North, y East, z Down (flat Earth)
  - Body: x out the nose, y out the right wing, z out the belly
  - Euler angles: 3-2-1 sequence (yaw psi, then pitch theta, then roll phi), radians
"""
import numpy as np

def quat_normalize(q):
    """Return q scaled to unit length. Call once per step, after RK4."""
    q = np.asarray(q, dtype=float)
    n = np.linalg.norm(q)
    return q / n

def dcm_from_quat(q):
    """
    Direction cosine matrix C_bn: rotates a vector from BODY axes to NED.
        v_ned  = C_bn   @ v_body
        v_body = C_bn.T @ v_ned
    Assumes q is unit length.
    """
    q0, q1, q2, q3 = q

    return np.array([
        [q0**2 + q1**2 - q2**2 - q3**2,  2*(q1*q2 - q0*q3),              2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3),              q0**2 - q1**2 + q2**2 - q3**2,  2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2),              2*(q2*q3 + q0*q1),              q0**2 - q1**2 - q2**2 + q3**2],
    ])

def euler_to_quat(phi, theta, psi):
    """
    3-2-1 Euler angles (radians) -> unit quaternion [q0, q1, q2, q3].
      phi   = roll  (about body x)
      theta = pitch (about body y)
      psi   = yaw   (about body z)
    """
    # every term uses HALF angles
    cphi, sphi = np.cos(phi / 2),   np.sin(phi / 2)
    cth,  sth  = np.cos(theta / 2), np.sin(theta / 2)
    cpsi, spsi = np.cos(psi / 2),   np.sin(psi / 2)

    q0 = cphi*cth*cpsi + sphi*sth*spsi
    q1 = sphi*cth*cpsi - cphi*sth*spsi
    q2 = cphi*sth*cpsi + sphi*cth*spsi
    q3 = cphi*cth*spsi - sphi*sth*cpsi
    return np.array([q0, q1, q2, q3])

    

def body_velocity(V, alpha, beta):
    """Inverse of air_angles: (V, alpha, beta) -> body velocity [u, v, w]."""
    u = V * np.cos(alpha) * np.cos(beta)
    v = V * np.sin(beta)
    w = V * np.sin(alpha) * np.cos(beta)
    return np.array([u, v, w])


def body_to_earth(vec, quat):
    """Rotate a vector from BODY axes to Earth (NED) axes."""
    return dcm_from_quat(quat) @ np.asarray(vec, dtype=float)


def earth_to_body(vec, quat):
    """Rotate a vector from Earth (NED) axes to BODY axes -- transpose of body_to_earth."""
    return dcm_from_quat(quat).T @ np.asarray(vec, dtype=float)


def quat_rate(quat, rate):
    """How fast the attitude quaternion changes, given body rates (q_dot = ½ q ⊗ [0, ω])."""
    q0, q1, q2, q3 = quat
    roll_rate, pitch_rate, yaw_rate = rate
    return 0.5 * np.array([
        -q1*roll_rate - q2*pitch_rate - q3*yaw_rate,
         q0*roll_rate - q3*pitch_rate + q2*yaw_rate,
         q3*roll_rate + q0*pitch_rate - q1*yaw_rate,
        -q2*roll_rate + q1*pitch_rate + q0*yaw_rate,
    ])