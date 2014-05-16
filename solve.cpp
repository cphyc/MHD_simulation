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
int main(){
  // Change the values of Ra and Pr
  Ra = 7e+3;
  Pr = 1; 
  dt = 0.50*pow(dz,2)/(4*max(1,Pr));

  // Create a new simulation
  Simulation s (20,0,0);

  // Init T with Tbound_top = 0, Tbount_bot = DT
  Temp T (T_init, 0, 0);

  // Init with null bound conditions
  Vort w (null_init, 0, 0);
  Stream phi (null_init, 0, 0);

  do {
      // Calculate the new T
      T.step();
      // Use this new T to calculate w
      w.step(T);
      // Use this new w to calculate phi
      phi.step(w);
  } while (s.iter(T, w, phi)); 
  return 0;
}
