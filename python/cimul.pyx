#!/usr/env python3
#cython: boundscheck=False
#cython: wraparound=False
import numpy as np
cimport numpy as np

from cpython cimport bool

cdef class Simulation

cdef class Vector:
    cdef np.ndarray G, G_old, field
    cdef str name
    cdef Simulation p
    def __init__(self, Simulation p):
        # save the pointer to the parent (dynamical)
        self.p = p
        # initial G = 0, G[k,n]
        self.G = np.zeros((self.p.Nz, self.p.NFourier), dtype="float64")
        # access via G[k][n]

    def step(self):
        cdet int end = self.p.Nz
        # save the old G
        self.G_old = self.G

        # compute the new one 
        self.compute_G()

        # new += dt/2*(3G-G_old)
        
        self.field[1:end-1,:] = (self.field[1:end-1,:] 
                              + self.p.dt/2*(3*self.G[1:end-1,:] - self.G_old[1:end-1,:])
                          ) 

        # conditions at top and bottom : null
        self.field[0 ,:] = 0
        self.field[end-1,:] = 0

    cpdef compute_G(self):
        pass

    def initial(self, str init_cond):
        if init_cond == 'null':
            self.field = np.zeros((self.p.Nz, self.p.NFourier))
        elif init_cond == "T":
            self.field = np.array([[T_0(n,k,self.p) for n in range(self.p.NFourier)]
                                   for k in range(self.p.Nz)])
        else:
            raise Exception("init_cond must be either `null` or `T`")

cdef class Temp(Vector):
    def __init__(self, parent):
        self.name = "T"
        super().__init__(parent)

    cpdef compute_G(self):
        cdet int end = self.p.Nz
        # compute G except for k = 0, Nz-1 and n = 0
        for n in range(1, self.p.NFourier):
            self.G[1:end-1,n] = ((self.field[:end-2,n]-2*self.field[1:end-1,n]+self.field[2:end,n]) 
                              * self.p.oodz2
                              - (n*self.p.pi/self.p.a)**2 
                              * self.field[1:end-1,n] )

cdef double sub_f (int k,int n, Simulation s):
    if k < (s.Nz-1):
        return s.oodz2
    else:
        return 1
cdef double dia (int k,int n, Simulation s):
    if 0 < k <(s.Nz-1):
        return (n*s.pi/s.a)**2 + 2*s.oodz2
    else:
        return 1

cdef double sup (int k,int n, Simulation s):
    if k > 0:
        return s.oodz2
    else:
        return 1
        
cdef class Vort(Vector):
    cdef np.ndarray wk1, wk2, sub
    def __init__(self, Simulation p):
        self.name = "ω"
        super().__init__(p)
        self.compute_wk()

    def compute_wk(self):
        cdet int end = self.p.Nz
        # init. the arrays:
        self.wk1 = np.zeros((self.p.Nz, self.p.NFourier))
        self.wk2 = np.zeros((self.p.Nz, self.p.NFourier))
        self.sub = np.zeros((self.p.Nz, self.p.NFourier))
        for n in range(1,self.p.NFourier):
            
            # tridiag. solver
            self.wk1[0,n] = 1/dia(0,n,self.p)
            self.wk2[0,n] = sup(0,n,self.p) * self.wk1[0,n]
            for k in range(1, self.p.Nz-1):
                self.wk1[k,n] = 1 /(dia(k,n,self.p)-sub_f(k,n,self.p)*self.wk2[k-1,n])
                self.wk2[k,n] = sup(k,n,self.p)*self.wk1[k,n]

            self.wk1[end-1,n] = 1/(dia(self.p.Nz-1,n,self.p)-sub_f(self.p.Nz-1,n,self.p)*self.wk2[end-2,n])

            self.sub[:,n] = [sub_f(k,n,self.p) for k in range(self.p.Nz)]

        
    cpdef compute_G(self):
        cdet int end = self.p.Nz
        cdef np.ndarray[np.float64_t, ndim=2] rhs
        rhs = self.p.psi.field
        # boundary conditions k=0, Nz-1 : psi = 0
        rhs[0, :] = 0
        rhs[end-1,:] = 0
        for n in range(1,self.p.NFourier):  
            # tridiag. solver
            self.field[0,n] = rhs[0,n]*self.wk1[0,n]
            for k in range(1, self.p.Nz):
                self.field[k,n] = (rhs[k,n] - self.sub[k,n]*self.field[k-1,n]*self.wk1[k,n])
                    
            for k in range(self.p.Nz-2, 0, -1):
                self.field[k,n] = self.field[k,n]-self.wk2[k,n]*self.field[k+1,n]
                    

            
cdef class Stream(Vector):
    def __init__(self, Simulation p):
        super().__init__(p)
        self.name = "ψ"

    cpdef compute_G(self):
        # compute G except for k=0, Nz-1 and n=0
        for n in range(1, self.p.NFourier):
            self.G[1:-1,n] = self.p.Pr*( self.p.Ra*n*self.p.pi/self.p.a*self.p.T.field[1:-1,n]
                                         + (self.field[:-2,n] - 2*self.field[1:-1,n] + self.field[2:,n])*self.p.oodz2
                                         - (n*self.p.pi/self.p.a)**2*self.field[1:-1,n] )

cdef class Simulation:
    cdef double dt_security, t
    cdef public double oodz2, pi, a, Re, Pr, Ra, dz, dt
    cdef public int Nz, NFourier
    cdef int freq_output, freq_critical_Ra, maxiter, niter
    cdef bool verbose
    cdef public Temp T
    cdef public Vort omega
    cdef public Stream psi
    cdef np.ndarray T_old, omega_old, psi_old

    def __init__(self, double Re= 1, double Pr= 1, double Ra= 1, double a = 1, int Nz= 100,
                 int NFourier= 50, double dt_security= 0.9,
                 int maxiter= 100, int freq_output= 10,
                 int freq_critical_Ra=50, bool verbose= False):
        self.Re = Re
        self.Pr = Pr
        self.Ra = Ra
        self.a = a
        self.Nz = Nz
        self.NFourier = NFourier
        self.dt_security = dt_security
        self.maxiter =  maxiter
        self.freq_output = freq_output
        self.freq_critical_Ra = freq_critical_Ra
        self.verbose = verbose

        # set the initial values
        self.t = 0
        self.niter = 0
        self.dz = 1./(1.0*self.Nz-1)

        # some usefull quantities
        self.oodz2 = 1/self.dz**2
        self.pi = np.pi

        # create the inner fields
        self.T = Temp(self)
        self.omega = Vort(self)
        self.psi = Stream(self)

        # previous fields for critical Ra number
        self.T_old = np.zeros((self.NFourier,1))
        self.omega_old = np.zeros((self.NFourier,1))
        self.psi_old = np.zeros((self.NFourier,1))

    def growth(self):
        ''' Calculate the log-growth rate and return a string containing
        all the growth rate'''
        cdef str out
        amp = lambda v: np.log(abs(v)) if v != 0 else 0
        gr = lambda new,old,n: str(amp(new.field[self.Nz//3,n])
                                   - amp(abs(old[n])))
        out = "".join([ gr(self.T, self.T_old,n) + "\t" +
                        gr(self.omega, self.omega_old,n) + "\t" +
                        gr(self.psi, self.psi_old,n) + "\t"
                        for n in range(self.NFourier) ])

        # save the arrays for next output
        self.T_old = self.T.field[self.Nz//3,:]
        self.omega_old = self.omega.field[self.Nz//3,:]
        self.psi_old = self.psi.field[self.Nz//3,:]

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
        
cdef double T_0 (int n,int k, s):
    if n > 0:
        return np.sin(s.pi*k*s.dz)
    else:
        return 1-k*s.dz
