from simul import *

if __name__ == '__main__':
    # create a new simulation
    s = Simulation(Re=5)

    # initial conditions psi(0) = 0, Omega(0) = 0
    s.psi.initial("null")
    s.omega.initial("null")
    # T_n(t=0) = sin(pi*k*dz) & T_0(t=0) = 1-k*dz
    s.T.initial(lambda n, k: T_0(n,k,s))

    # main loop over time
    while s.step():
        s.T.step()
        s.psi.step()
        s.omega.step()

    del s
        
