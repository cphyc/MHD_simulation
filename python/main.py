#!/usr/env python3
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle
import ipdb

class Vector(object):
    def __init__(self, parent):
        # save the pointer to the parent (dynamical)
        self.p = parent
        # initial G = 0, G[k,n]
        self.G = np.zeros((self.p.Nz, self.p.NFourier), dtype="float64")
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
        if init_cond == 'null':
            self.field = np.zeros((self.p.Nz, self.p.NFourier))
        else: 
            try:
                self.field = np.array([[init_cond(n,k) for n in range(self.p.NFourier)]
                                       for k in range(self.p.Nz)])
            except TypeError:
                raise Exception("init_cond must be either `null` or a 2-arguments callable")

class Temp(Vector):
    name = "T"
    def compute_G(self):
        # compute G except for k = 0, Nz-1 and n = 0
        for n in range(1, self.p.NFourier):
            self.G[1:-1,n] = ((self.field[:-2,n]-2*self.field[1:-1,n]+self.field[2:,n]) 
                              * self.p.oodz2
                              - (n*self.p.pi/self.p.a)**2 
                              * self.field[1:-1,n] )
        
class Vort(Vector):
    name = "ω"
    def __init__(self, parent):
        super().__init__(parent)
        self.compute_wk()

    def compute_wk(self):
        # init. the arrays:
        self.wk1 = np.zeros((self.p.Nz, self.p.NFourier))
        self.wk2 = np.zeros((self.p.Nz, self.p.NFourier))
        self.sub = np.zeros((self.p.Nz, self.p.NFourier))
        for n in range(1,self.p.NFourier):
            # save some usefull functions
            sub_f = lambda k : -self.p.oodz2 if k<self.p.Nz-1 else 1
            dia = lambda k : (n*self.p.pi/self.p.a)**2 + 2*self.p.oodz2 if 0<k<self.p.Nz-1 else 1
            sup = lambda k : -self.p.oodz2 if k>0 else 1
            
            # tridiag. solver
            self.wk1[0,n] = 1/dia(0)
            self.wk2[0,n] = sup(0) * self.wk1[0,n]
            for k in range(1, self.p.Nz-1):
                self.wk1[k,n] = 1 /(dia(k)-sub_f(k)*self.wk2[k-1,n])
                self.wk2[k,n] = sup(k)*self.wk1[k,n]

            self.wk1[-1,n] = 1/(dia(self.p.Nz-1)-sub_f(self.p.Nz-1)*self.wk2[n,-2])

            self.sub[:,n] = [sub_f(k) for k in range(self.p.Nz)]

        
    def compute_G(self):
        rhs = self.p.psi.field
        # boundary conditions k=0, Nz-1 : psi = 0
        rhs[0, :] = 0
        rhs[-1,:] = 0
        for n in range(1,self.p.NFourier):  

            # tridiag. solver
            self.field[0,n] = rhs[0,n]*self.wk1[0,n]
            for k in range(1, self.p.Nz):
                self.field[k,n] = (rhs[k,n] - self.sub[k,n]*self.field[k-1,n]*self.wk1[k,n])
                    
            for k in range(self.p.Nz-2, 0, -1):
                self.field[k,n] = self.field[k,n]-self.wk2[k,n]*self.field[k+1,n]
                    

            
class Stream(Vector):
    name = "ψ"
    def compute_G(self):
        # compute G except for k=0, Nz-1 and n=0
        for n in range(1, self.p.NFourier):
            self.G[1:-1,n] = self.p.Pr*( self.p.Ra*n*self.p.pi/self.p.a*self.p.T.field[1:-1,n]
                                         + (self.field[:-2,n] - 2*self.field[1:-1,n] + self.field[2:,n])*self.p.oodz2
                                         - (n*self.p.pi/self.p.a)**2*self.field[1:-1,n] )

class Simulation(object):
    param_list = {'Re': 1, 'Pr': 1, 'Ra': 1, 'a' : 1, 'Nz': 100,
                  'NFourier': 50, 'dt_security': 0.9,
                  'maxiter': 100, 'freq_output': 10,
                  'output_file': "dump.log"}

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

        # set the initial values
        self.t = 0
        self.niter = 0
        self.dz = 1/(self.Nz-1)

        # some usefull quantities
        self.oodz2 = 1/self.dz**2
        self.pi = np.pi

        # create the inner fields
        self.T = Temp(self)
        self.omega = Vort(self)
        self.psi = Stream(self)

        # T previous for critical Ra number

    def __del__(self):
        pass

    def step(self):
        # eventually output
        if self.niter % self.freq_output == 0:
            self.dump()

        # eventually calculate the d-ln term for the critical Ra
        if self.niter % self.freq_cirtical_Ra == 0:
            self.T_old = self.T.field[self.Nz//3,n]
            output = "# T growth : \t"
            for n in range(self.NFourier):
                output += ""
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

    def dump(self):
        output = "#k\t"
        for n in range(self.NFourier):
            o = "{T}_{n}\t{w}_{n}\t{psi}_{n}\t".format(T=self.T.name,
                                                       w=self.omega.name,
                                                       psi=self.psi.name,
                                                       n=n)
            output += o
        output += "\n"
        for k in range(self.Nz):
            output += str(k) + "\t"
            for n in range(self.NFourier):
                l = "{T}\t{w}\t{psi}\t".format(T=self.T.field[k,n],
                                               w=self.omega.field[k,n],
                                               psi=self.psi.field[k,n])
                output += l
            output += "\n"
        print(output)

    def CFL(self):
        # dt < (dz)^2/4 or (dz)^2/(4Pr) if Pr > 1
        self.dt = self.dt_security * self.dz**2/(4*max(1,self.Pr))
        
def T_0 (n,k,s):
    if n > 0:
        return np.sin(s.pi*k*s.dz)
    else:
        return 1-k*s.dz

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

    #del s
        
        
