#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <string>
#include "classes.hpp"
#include "misc.hpp"

Simulation::Simulation (int save_freq, int save_freq_growth, 
			int max_niter, double max_time,
			double init_time, int init_niter,
			double dt_security, double dt_additional_cfl_security,
			double init_dt){
  this->save_freq = save_freq;
  this->t = init_time; 
  this->niter = init_niter;
  this->save_freq_growth = save_freq_growth;
  this->max_niter = max_niter;
  this->max_time = max_time;
  this->dt_security = dt_security;
  this->dt_additional_cfl_security = dt_additional_cfl_security;
  this->dt = init_dt;
  T = new Temp(0, 0, this);
  w = new Vort(0, 0, this);
  psi = new Stream(0, 0, this);
}


// Get an array of T, w, psi and return a string containing :
// T0 w0 psi0 T1 w1 psi1 ... Tn wn psin
void Simulation::dump(){
  // make a FFT transform
  do_FFT();
  
  // write the header
  std::cout << "# time\t" << t << std::endl;
  std::cout << "# niter\t" << niter << std::endl;

  // count the columns
  std::cout <<  "# Cols ";
  for (int nc = 1; nc < (NMAX-1)*4+1; nc++)
    std::cout << nc << "\t";
  std::cout << std::endl;

  // print the title
  std::cout << "# Header k" ;
  for (int kx = 0; kx < NX; kx++)
    std::cout << "\tT(x=" << kx << ")\tω(x=" << kx << ")\tψ(x=" << kx << ")" ;
  
  // get back to the next line
  std::cout << "\n";

  // write the data (one line = one k)
  for (int k = 0; k < NZ; k++){
    std::cout << k << "\t";
    
    for (int kx = 0; kx < NX; kx++){
      //std::cout << "#k:n :\t" << k <<":"<< n << "\n";
      std::cout << (T->real_val[kx][k]) << "\t" ;
      std::cout << (w->real_val[kx][k]) << "\t";
      std::cout << (psi->real_val[kx][k]) << "\t";
    }
    // ended one n, we get back to the next line for the next k
    std::cout<< "\n";
  }
  std::cout << "\n";
}

void Simulation::cfl() {
  // Check for the cfl each cfl_freq
  if (cfl_freq == 0 || niter % cfl_freq == 0) {
    do_FFT();

    dt_old = dt;
    bool is_along_x;
    double vmax = 0, vx = 0, vz = 0, new_dt = 0;
    // get the maximum of the speed (derivative of psi)
    for (int z = 1; z < NX-1; z++){
      for (int x = 1; x < NX-1; x++){
	vx = (this->psi->real_val[x+1][z] - this->psi->real_val[x-1][z])/(2*DX);
	vz = (this->psi->real_val[x][z+1] - this->psi->real_val[x][z-1])/(2*DZ);
	if (vx > vmax && vx >= vz) {
	  vmax = vx;
	  is_along_x = true;
	}
	else if (vz > vmax && vz > vx) {
	  vmax = vz;
	  is_along_x = false;
	}
      }
    }
    //dt = dt_security * pow(DZ,2)/4*max(1,Pr);
    if (is_along_x)
      new_dt= dt_additional_cfl_security * dt_security * DX/vmax;
    else
      new_dt= dt_additional_cfl_security * dt_security * DZ/vmax;

    double dt_non_cfl = dt_security*(DZ*DZ)/(4*max(Pr,1));
    if (dt_non_cfl < new_dt)
      this->dt= dt_non_cfl;

    std::cout << "# CFL: Updated dt from " << dt_old << " to " << dt << std::endl;
  }

}

void Simulation::unset_FFT() {
  w->unset_FFT();
  psi->unset_FFT();
  T->unset_FFT();
}

void Simulation::do_FFT() {
  w->FFT();
  psi->FFT();
  T->FFT();
}
// Iterator of the simulation. Stops when the MAX_TIME or MAX_ITER is reached.
bool Simulation::iter(){
  // Compute dt
  cfl();

  // Say that the FFT has not been done
  unset_FFT();

  // Save every save_freq except if <= 0
  if (save_freq > 0 && niter % save_freq == 0) {
    this->dump();
  }
  // Save the growth rate every save_freq_growth except if <=0

  if (save_freq_growth > 0 && niter % save_freq_growth == 0){
    // Growth rate calculated at a third of the heigth
    int third = NZ / 3;

    // Save the initial fields
    if (niter == 0) {
      for (int n = 1; n < NMAX; n++){
	old_T[n] = T->val[n][third];
	old_w[n] = w->val[n][third];
	old_psi[n] = psi->val[n][third];
      }
    }
	  
    // Iterate over all the modes
    for (int n = 1; n < NMAX; n++){
      // compute log | T_old | / | T |
      double gr_T = log( abs(T->val[n][third])) - log(abs(old_T[n]) );
      double gr_w = log(abs(w->val[n][third])) - log(abs(old_w[n]));
      double gr_psi = log(abs(psi->val[n][third])) - log(abs(old_psi[n]));

      // output it
      std::cout << "# Growth rate :\t" << n << "\t";
      std::cout << niter << "\t" << gr_T;
      std::cout << "\t" << gr_w << "\t" << gr_psi << "\n";

      // save it for the next output
      old_T[n] = T->val[n][third];
      old_w[n] = w->val[n][third];
      old_psi[n] = psi->val[n][third];      
    }
  }
  
  // move of one time step and check we can continue
  niter ++;
  t += dt;
  if ((max_time > 0 && t > max_time) || (max_niter > 0 && niter > max_niter)) return false;
  return true;
}
