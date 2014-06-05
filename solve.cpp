#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <string>
#include "classes.hpp"
#include "misc.hpp"
 
using namespace std;

int main(){
  // Create a new simulation
  // arguments are : save_freq, save_freq_growth, 
  //		     max_niter, max_time,
  //		     init_time, init_niter, dt_security
  // all optionnalsg

  // The simulation creates T,w,psi
  Simulation s (1);
  // Time and output control
  s.cfl_freq = 10;
  s.save_freq = 500;
  s.max_niter = 100000;
  s.dt = 1e-8;
  s.dt_old = s.dt;
  s.dt_security = 0.8;
  s.dt_additional_cfl_security = 0.5;

  // Simulation parameter
  s.Ra = 1e6;
  s.Re = 1;
  s.Pr = 0.5;

  for (int n = 0; n < NMAX; n++) {
    for (int k = 0; k < NZ; k++) {
      s.T->val[n][k] = sin(PI*k*DZ);
    }
  }

  // Iterate at least once over time
  do {
    // Compute the difference
    s.psi->compute();
    s.w->compute();
    s.T->compute();

    // Add it
    s.psi->add();
    s.w->add();
    s.T->add();

    // Boundary conditions
    s.psi->boundaries();
    s.w->boundaries();
    s.T->boundaries();

    // Compute FFT
    s.psi->FFT();
    s.w->FFT();
    s.T->FFT();
  } while (s.iter());

  return 0;
}
