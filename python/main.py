#!/usr/env python3
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle
        

class Vector(object):
    def __init__(self, parent):
        # save the pointer to the parent (dynamical)
        self.p = parent
        # initial G = 0, G[k,n]
        self.G = np.zeros((parent.Nz, parent.NFourier), dtype="float64")
        # access via G[k][n]

    def step(self):
        # save the old G
        self.G_old = self.G

        # compute the new one 
        self.compute_G()

        # new += dt/2*(3G-G_old)
        self.field = self.field + self.p.dt/2*(3*self.G - self.G_old) 

    def compute_G(self):
        raise Exception("Vector class is a base class, not supposed to be "+
                        "used like that")

    def initial(self, init_cond):        
        # set the initial condition
        self.field = np.fromfunction(init_cond, (self.p.Nz, self.p.NFourier))

class Temp(Vector):
    def compute_G(self):
        # compute G except for k = 0, Nz-1 and n = 0
        for n in range(1, self.p.NFourier):
            self.G[1:-1,n] = ((self.field[:-2,n]-2*self.field[1:-1,n]+self.field[2:,n]) 
                              * self.p.oodz2
                              - (n*self.p.pi/self.p.a)**2 
                              * self.field[1:-1,n] )
        
class Vort(Vector):
    def compute_G(self):
        for n in range(1,self.p.NFourier):
            # save some usefull functions
            sub = lambda k : -self.p.oodz2 if k<Nz-1 else 1
            dia = lambda k : (n*self.p.pi/self.p.a)**2 + 2*self.p.oodz2 if 0<k<Nz-1 else 1
            sup = lambda k : -self.p.oodz2 if k>0 else 1
            self.rhs = self.p.omega.field[:,n]
        
            # boundary conditions k=0, Nz-1 : phi = 0
            for n in range(self.p.NFourier):
                self.rhs[ 0] = 0
                self.rhs[-1] = 0

            # tridiag. solver
            wk1[0] = 1/dia(0)
            wk2[0] = sup(0) * wk1[0]
            for k in range(1,self.p.Nz-1):
                wk1[k] = 1 /(dia(k)-sub(k)*wk2(k-1))
                wk2[k] = sup(k)*wk1[k]
            wk1[-1] = 1/(dia(self.p.Nz-1)-sub(self.p.nz-1)*wk2(self.p.nz-2))

            sol[0] = rhs[0]*wk1(0)
            for k in range(1, self.p.Nz-1,  1):
                sol[k] = (rhs(k) - sub(k)*sol[k-1]*wk1(k))
            for k in range(self.p.Nz-2, 0, -1):
                sol[k] = sol[k]-wk2(k)*sol[k+1]

            self.field[n] = sol

            
class Stream(Vector):
    pass

class Simulation(object):
    param_list = {'Re': 1, 'Pr': 1, 'a' : 1, 'Nz': 100,
                  'NFourier': 50, 'dt_security': 0.9,
                  'maxiter': 10000}

    def __init__(self, *args, **kargs):
        # save the default parameters
        for param, value in self.param_list.items():
            setattr(self, param, value)

        # override if necessary
        for param, value in kargs.items():
            if param not in self.param_list:
                raise Exception("`%s' not recognized" % param)
            else:
                setattr(self, param, value)

        # create the inner fields
        self.T = Temp(self)
        self.omega = Vort(self)
        self.phi = Stream(self)

        # set the initial values
        self.t = 0
        self.niter = 0
        self.dz = 1/(self.Nz-1)

        # some usefull quantities
        self.oodz2 = 1/self.dz**2
        self.pi = np.pi

    def step(self):
        # get the max timestep
        self.CFL()

        # increase the time, the iteration
        self.t += self.dt
        self.niter += 1

        # check that the end is not reached
        if self.niter > self.maxiter:
            return False
        else:
            return True

    def CFL(self):
        # dt < (dz)^2/4 or (dz)^2/(4Pr) if Pr > 1
        self.dt = self.dt_security * self.dz**2/(4*max(1,self.Pr))
        
if __name__ == '__main__':
    # create a new simulation with border conditions
    s = Simulation(Re=5)

    # initial conditions T = np.sin(pi*z(k)) with z(k) = k*dz (k in [0,Nz-1])
    # and T0(k) = 1 - k*dz
    s.phi.initial(lambda _n,_k: 0)
    s.omega.initial(lambda _n,_k: 0)
    s.T.initial(lambda n, k: (n!=0) * np.sin(s.pi*k*s.dz) + n==0*(1-k*s.dz))

    # main loop over time
    while s.step():
        s.T.step()
        s.omega.step()
        s.phi.step()

    #     s.dump()
        
        
