#!/usr/env python3
try:
    import numpypy as np
except:
    import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle
import ipdb

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    nf = len(a) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy the array
    for it in xrange(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]
 
    xc = ac
    xc[-1] = dc[-1]/bc[-1]
 
    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
 
    del bc, cc, dc # delete variables from memory
 
    return xc

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
        self.field[1:-1] = (self.field[1:-1] 
                            + self.p.dt/2*(3*self.G[1:-1] - self.G_old[1:-1])
                        ) 

        # conditions at top and bottom : null
        self.field[0 ,:] = 0
        self.field[-1,:] = 0

    def compute_G(self):
        raise Exception("Vector class is a base class, not supposed to be "+
                        "used like that")

    def initial(self, init_cond):
        if init_cond == 'null':
            self.field = np.zeros((self.p.Nz, self.p.NFourier))
        elif init_cond == "T":
            self.field = np.array([[T_0(n,k,self.p) for n in range(self.p.NFourier)]
                                   for k in range(self.p.Nz)])
        else:
            raise Exception("init_cond must be either `null` or `T`")

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

            self.wk1[-1,n] = 1/(dia(self.p.Nz-1)-sub_f(self.p.Nz-1)*self.wk2[-2,n])

            self.sub[:,n] = [sub_f(k) for k in range(self.p.Nz)]

    def step(self):
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
                  'freq_critical_Ra':50, 'verbose': False}

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

        # previous fields for critical Ra number
        self.T_old = np.zeros((self.NFourier,))
        self.omega_old = np.zeros((self.NFourier,))
        self.psi_old = np.zeros((self.NFourier,))

    def __del__(self):
        pass

    def growth(self):
        ''' Calculate the log-growth rate and return a string containing
        all the growth rate'''
        
        amp = lambda v: np.log(abs(v)) if v != 0 else 0
        gr = lambda new,old,n: str(amp(new.field[self.Nz//3,n])
                                   - amp(abs(old[n])))
        out = "".join([ gr(self.T, self.T_old,n) + "\t" +
                        gr(self.omega, self.omega_old,n) + "\t" +
                        gr(self.psi, self.psi_old,n) + "\t"
                        for n in range(self.NFourier) ])

        # save the arrays for next output
        ipdb.set_trace()
        self.T_old = self.T.field[self.Nz//3,:].copy()
        self.omega_old = self.omega.field[self.Nz//3,:].copy()
        self.psi_old = self.psi.field[self.Nz//3,:].copy()

        return out+"\n"

    def step(self):
        # eventually output
        if self.verbose and self.niter % self.freq_output == 0:
            self.dump()

        # eventually calculate the d-ln term for the critical Ra
        if self.verbose and self.niter % self.freq_critical_Ra == 0 :
            output = "# growth : \t"
            output+= "".join([
                "{T.name}_{n}\t{w.name}_{n}\t{psi.name}_{n}\t".format(T=self.T,
                                                                    w=self.omega,
                                                                    psi=self.psi,
                                                                    n=n)
                for n in range(self.NFourier)])
            output+= "\n"
            output+= "# growth : \t"
            output+= self.growth()

            print(output)

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

    del s
