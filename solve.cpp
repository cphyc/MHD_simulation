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
  Simulation s (20, 0, 1000);

  s.Ra = 1;
  s.Re = 1;
  s.Pr = 1;

  // Iterate at least once over time
  do {    
    // Other way : calculate linear first psi, then w, then T
    s.psi->step();
    s.w->step();
    s.T->step();

  } while (s.iter()); 

  return 0;
}
