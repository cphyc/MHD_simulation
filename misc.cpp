#include <cmath>
#include <iostream>
#include "misc.hpp"

/* Solve a matrix equation of the form :
   M*A = B with M tridiagonal and M,B known */
void triDiSolve(vector &rhs, vector &sol, vector &sub,
		vector &dia, vector &sup){
  vector wk1, wk2;

  wk1[0] = 1.0/dia[0];
  wk2[0] = sup[0]*wk1[0];

  for(int i=1; i<NZ-1; i++){
    wk1[i] = 1.0 / (dia[i]-sub[i]*wk2[i-1]);
    wk2[i] = sup[i]*wk1[i];
  }

  wk1[NZ-1] = 1.0 / (dia[NZ-1]-sub[NZ-1]*wk2[NZ-2]);

  sol[0] = rhs[0]*wk1[0];

  for(int i=1; i<NZ-1; i++){
    sol[i] = (rhs[i]-sub[i]*sol[i-1]) * wk1[i];
  }
  for(int i=NZ-2; i>0; i--){
    sol[i] = sol[i] - wk2[i]*sol[i+1];
  }
}

// // Vector full of zeros
// vector null_vector() {
//   vector n;
//   for (int i = 0; i<NZ; i++){
//     n[i] = 0;
//   }
//   return n;
// }


// Initialization functions

double T_init (int n, int k) { return sin(PI*k*DZ);}
double T_0 (int k) { return DT*(1-DZ*k); }

double null_init(int n, int k) {return 0;}
double null_0 (int k) { return 0; }

double abs(double x) {
  if (x >= 0) return x;
  else return -x;
}

void pp(double val) { std::cout << val << "\t"; }
