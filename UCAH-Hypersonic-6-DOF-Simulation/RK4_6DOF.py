def rk4(f, t, state, dt):
 
    k1, outputs = f(t,          state)
    k2, _       = f(t + dt/2,   state + dt/2 * k1)
    k3, _       = f(t + dt/2,   state + dt/2 * k2)
    k4, _       = f(t + dt,     state + dt   * k3)

    return state + dt/6 * (k1 + 2*k2 + 2*k3 + k4), outputs
