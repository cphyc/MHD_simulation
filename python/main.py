from simul import *
import sys

if __name__ == '__main__':
    # create a new simulation
    print("# growth",sys.argv[1])
    s = Simulation(verbose=True, Ra=float(sys.argv[1]), maxiter=10000,
                   NFourier=10, Nz = 100)

    # initial conditions psi(0) = 0, Omega(0) = 0
    s.psi.initial("null")
    s.omega.initial("null")
    # T_n(t=0) = sin(pi*k*dz) & T_0(t=0) = 1-k*dz
    s.T.initial("T")

    # main loop over time
    while s.step():
        s.T.step()
        s.omega.step()
        s.psi.step()

    # remove if debugging !
    del s
