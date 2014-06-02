#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <string>
#include "classes.hpp"
#include "misc.hpp"

using namespace std;

// Initialize Ra and Pr
double Ra;
double Pr;
double dt;
int MAX_NITER = 10000;

int main(){
  // Change the values of Ra and Pr
  Ra = 1000;
  Pr = 1; 
  dt = 0.50*pow(dz,2)/(4*max(1,Pr));

  // Create a new simulation
  Simulation s (20,0,0);

  // Init T with Tbound_top = 0, Tbount_bot = DT
  Temp T (T_init, 0, 0);

  // Init with null boundary conditions
  Vort w (null_init, 0, 0);
  Stream phi (null_init, 0, 0);

  // Iterate at least once over time
  do {    
    // Other way : calculate first phi, then w, then T
    phi.step(w);
    w.step(T);
    T.step();

  } while (s.iter(T, w, phi)); 

  return 0;
}
