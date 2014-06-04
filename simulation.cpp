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
			double dt_security){
  this->save_freq = save_freq ;
  this->t = init_time; 
  this->niter = init_niter;
  this->save_freq_growth = save_freq_growth;
  this->max_niter = max_niter;
  this->max_time = max_time;
  this->dt_security = dt_security;
  T = new Temp(T_init, 0, 0, this);
  w = new Vort(null_init, 0, 0, this);
  psi = new Stream(null_init, 0, 0, this);
}


// Get an array of T, w, psi and return a string containing :
// T0 w0 psi0 T1 w1 psi1 ... Tn wn psin
void Simulation::dump(){
  // write the header
  std::cout << "\n\n#Header t:\t" << t << "\t" << dt_security*pow(DZ,2)/(4*max(1,Pr));

  // count the columns
  std::cout << std::endl << "#";
  for (int n = 1; n < NMAX; n++)
    std::cout << n << "\t";
  std::cout << std::endl;

  // print the title
  std::cout << "k" ;
  for (int n = 1; n < NMAX; n++)
    std::cout << "\tT_" << n << "\tw_" << n << "\tpsi_" << n ;
  
  // get back to the next line
  std::cout << "\n";

  // write the data (one line = one k)
  for (int k = 0; k < NZ; k++){
    std::cout << k << "\t";
    
    for (int n = 0; n < NMAX; n++){
      //std::cout << "#k:n :\t" << k <<":"<< n << "\n";
      std::cout << (T->val[n][k]) << "\t" ;
      std::cout << (w->val[n][k]) << "\t";
      std::cout << (psi->val[n][k]) << "\t";
    }
    // ended one n, we get back to the next line for the next k
    std::cout<< "\n";
  }
}

void Simulation::cfl() {
  dt = dt_security * pow(DZ,2)/4*max(1,Pr);
}

// Iterator of the simulation. Stops when the MAX_TIME or MAX_ITER is reached.
bool Simulation::iter(){
  // Compute dt
  cfl();
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
