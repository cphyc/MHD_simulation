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
  // all optionnals

  // The simulation creates T,w,psi
  Simulation s (1);
  s.max_niter = 1000;
  s.Ra = 1;
  s.Re = 1;
  s.Pr = 1;

  // Iterate at least once over time
  while (s.iter()) {    
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

  }

  return 0;
}
