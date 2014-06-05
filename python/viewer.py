#!/usr/bin/python3
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

class Snapshot(object):
    def to_real_grid(self, NZ, a):
        ''' Converts the grid to the real space.'''
        NX = int(NZ*a)
        # The arrays are arranged as array[x,z]
        # Create the array of the same size as NX
        T_coeffs, w_coeffs, psi_coeffs = [np.zeros((NX,NZ))]*3
        T,w,psi = [np.zeros((NX,NZ))]*3			
    
        # Get the first nmode coeffs of this array and fill the rest with 0
        T_coeffs[:self.nmode] = self.data[1::3]
        w_coeffs[:self.nmode] = self.data[2::3]
        psi_coeffs[:self.nmode] = self.data[3::3]

		
        for k in range(self.nz):
            T[:,k] = np.real(np.fft.ifft(T_coeffs[:,k]))
            w[:,k] = np.real(np.fft.ifft(w_coeffs[:,k]))
            psi[:,k] = np.real(np.fft.ifft(psi_coeffs[:,k]))
		

        return (T, w, psi)

def to_object(h,d):
    out = Snapshot()
    # In the header : one info per line
    for line in h:
        if "niter" in line or "time" in line:
            setattr(out, line[1], line[2])
        else:
            setattr(out, line[1], line[2:])

    setattr(out, "header", h)
    setattr(out, "data", d)
    # save the columns individually
    for column, title in zip(out.data, out.Header):
        setattr(out, title, column)

    # save the arrays for each field
    setattr(out, "T", out.data[1::3])
    setattr(out, "w", out.data[2::3])
    setattr(out, "psi", out.data[3::3])
    return out

def chunkify(path = "../output.log"):
    '''Opens the path and cut it where it is supposed to be cut, then slice it towards
    a more readable format and returns an object that has the magic attribute "to_real_grid"
    that returns all the fields in their real form (not spectral in x-direction).'''
    f = open(path, "r")
    data = []
    infos = []
    cols = []
    header = []
    # Iterate over the lines until we reach two empty lines
    for l in f.readlines():
        line = l.split()
        
        if len(line) == 0:
            yield to_object(np.array(header),
                            np.array(data).transpose().astype(float))
            data = []
            header = []
        elif "#" in line:
            header.append(line)
        else:
            data.append(line)

if __name__ == "__main__":
    import sys

    # for snap in chunkify(sys.argv[1]):
    #     T, _,_ = snap.to_real_grid(100,1)
    #     niter = snap.niter
    #     time = snap.time
    #     plt.imshow(T)
    #     print("{}/{:0>4}.png".format(sys.argv[2], niter) )
    #     plt.title("Time $t=%f$, $n_{iter} = %d$" % (time, niter))

    #     plt.savefig("{}/{:0^4}.png".format(sys.argv[2], niter) )
    #     plt.cla()
        
    def save_snap(snap):
        niter = int(snap.niter)
        time = float(snap.time)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(snap.T)
        cbar = fig.colorbar(cax)

        print("{}/{:0>4}.png".format(sys.argv[2], niter) )
        plt.title("Time $t=%f$, $n_{iter} = %d$" % (time, niter))

        fig.savefig("{}/{:0>4}.png".format(sys.argv[2], niter) )

        fig.delaxes(fig.axes[1])

    # real_snap = ( (snap.to_real_grid(100,3), snap.niter, snap.time) for snap in chunkify(sys.argv[1]) )
    # try:
    try:
        nprocess = int(sys.argv[3])
    except:
        nprocess = 0

    if nprocess > 0:
        from multiprocessing import Pool
        print ("Multiprocessing with %d cores" % nprocess)
        pool = Pool(processes=nprocess)
        data = pool.map(save_snap, chunkify(sys.argv[1]))
    else:
         print ("No multiprocessing :(")
         for snap in chunkify(sys.argv[1]):
             save_snap(snap)

    import os
    os.system("mencoder 'mf://%s/*.png' -mf type=png:fps=10" % sys.argv[2]+
              " -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s/animation.mpg" % sys.argv[2])
