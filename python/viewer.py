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
	pass

def to_object(h,d):
	out = Snapshot()
	setattr(out, "time", float(h[0][2]))
	setattr(out, "niter", int(h[0][4]))
	setattr(out, "header", h)
	setattr(out, "data", d)
    # save the columns individually
	for column, title in zip(d,h[2][1:]):
		setattr(out, title, column)
		
	# save the number of modes
	setattr(out, "nmode", len(h[2][2:])//3)
	setattr(out, "nz", out.data.shape[1])

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
        
		if line == []:
			if len(data) > 0:
				yield to_object(np.array(data[:3]),
								np.array(data[3:]).transpose().astype(float))
				data = []
			else:
				chunkify()
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
        T,w,psi = snap.to_real_grid(100,1)
        niter = snap.niter
        time = snap.time
            
        plt.imshow(T)
        print("{}/{:0>4}.png".format(sys.argv[2], niter) )
        plt.title("Time $t=%f$, $n_{iter} = %d$" % (time, niter))

        plt.savefig("{}/{:0>4}.png".format(sys.argv[2], niter) )
        plt.cla()

    # real_snap = ( (snap.to_real_grid(100,3), snap.niter, snap.time) for snap in chunkify(sys.argv[1]) )
    # try:
    from multiprocessing import Pool
    nprocess = int(sys.argv[3])
    print ("Multiprocessing with %d cores" % nprocess)
    pool = Pool(processes=nprocess)
    data = pool.map(save_snap, chunkify(sys.argv[1]))
    # except:
    #     print ("No multiprocessing :(")
    #     for snap in real_snap:
    #         save_snap(snap)

    import os
    os.system("mencoder 'mf://%s/*.png' -mf type=png:fps=10" % sys.argv[2]+
              " -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s/animation.mpg" % sys.argv[2])
