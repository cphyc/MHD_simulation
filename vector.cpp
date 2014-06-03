#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <string>
#include "constant.cpp"
#include "classes.hpp"
#include "misc.hpp"

// Takes a function as argument and fill the initial value 
// thanks to the function
Vector::Vector(double (*init)(int, int), double top, double bot){
  for (int n = 1; n < NMAX; n++){
    for (int k = 0; k < NZ ; k++){
      val[n][k] = init(n, k);
    }
  }
  this->top_boundary = top;
  this->bot_boundary = bot;
}

// set boundary conditions at k = NZ-1 and k = 0
void Vector::boundaries(double (*val)(int)){
  // set boundaries
  for (int n = 1; n < NMAX; n++){
    this->val[n][NZ-1] = this->top_boundary;
    this->val[n][0   ] = this->bot_boundary;
  }

  // set mode 0
  for (int k = 0; k < NZ ; k ++)
    this->val[0][k] = val(k);
}

 
void Temp::compute_G(int n, int k){
  G_old[n][k] = G[n][k];

  G[n][k] = (val[n][k+1]
	     - 2*val[n][k]+val[n][k-1]) * oodz_2 
    - pow(n * PI / a, 2) * val[n][k];
}

// Iterates over all modes and compute the new linear_step
void Temp::linear_step(Vector psi){
  double nlt;
  for (int n = 0; n < NMAX; n++) {
    for (int k = 1; k < NZ-1; k++){
      // compute G
      this->compute_G(n,k);
	
      // Non linear term
      nlt = this->non_linear_term(psi, n, k);

      // New value
      this->val[0][k] += nlt + dt/2 * (3*G[n][k] - G_old[n][k]);
    }
  }
  // set the boundary conditions and apply T_0 for 0th mode
  this->boundaries(T_0);
}

double Temp::non_linear_term(Vector psi, int n, int k){
  // nl term from T0
  double nlt = 0;
  if (n == 0) {
    for (int n2 = 1; n2 < NMAX; n2++) {
      nlt += -PI/(2*a)*n2*( ( psi.val[n2][k+1]-psi.val[n2][k-1])/(2*dz)*this->val[n2][k]
			    + psi.val[n2][k]*(this->val[n2][k+1] - this->val[n2][k-1]/(2*dz)));
    }
    return nlt;
  }
  else {
    int n3;

    // Contribution for n' = 0, n'' = n
    nlt = -n*PI/a * psi.val[n][k]*(this->val[0][k+1] - this->val[0][k-1])/(2*dz);

    // nl term from other n', n''
    for (int n2 = 1; n2 < NMAX; n2++) {
      n3 = n - n2;
      if (n3 > 0 && n3 < NMAX) {
	nlt += -n*PI/a*psi.val[n][k]*(this->val[0][k+1] - this->val[0][k-1])
	  - PI/(2*a) * (-n2*(psi.val[n3][k+1] - psi.val[n3][k-1])/(2*dz) * this->val[n2][k]
			+n3*psi.val[n3][k] * (this->val[n2][k+1] - this->val[n2][k-1])/(2*dz));
      }
      n3 = n2 -n;
      if (n3 > 0 && n3 < NMAX) {
	nlt += -PI/(2*a)*(n2*(psi.val[n2][k+1] - psi.val[n3][k-1])/(2*dz) * this->val[n2][k]
			  +n3*psi.val[n3][k] * (this->val[n2][k+1] - this->val[n2][k-1])/(2*dz));
      }
      n3 = n+n2;
      if (n3 > 0 && n3 < NMAX) {
	nlt += -PI/(2*a)*(n2*(psi.val[n2][k+1] - psi.val[n3][k-1])/(2*dz) * this->val[n2][k]
			  +n3*psi.val[n3][k] * (this->val[n2][k+1] - this->val[n2][k-1])/(2*dz));
      }
    }
    return nlt;
  }
}

void Vort::compute_G(Vector T, int n, int k){

  G_old[n] = G[n];
  G[n][k] = Ra*Pr*(n * PI / NX) * T.val[n][k] 
    + Pr*((val[n][k+1] - 2*val[n][k]+val[n][k-1]) * oodz_2
	  - pow(n * PI / NX, 2) * val[n][k]);
}

  
void Vort::linear_step(Vector T, Vector psi){
  double nlt;
  for (int n = 0; n < NMAX; n++) {
    for (int k = 1; k < NZ-1; k++){
      // compute G
      this->compute_G(T, n, k);
	
      // Non linear term
      nlt = this->non_linear_term(psi, n, k);

      // New value
      this->val[0][k] += nlt + dt/2 * (3*G[n][k] - G_old[n][k]);
    }
  }
  // set the boundary conditions and apply 0 for 0th mode
  this->boundaries(null_0);
}

double Vort::non_linear_term(Vector w, int n, int k){
  if (n == 0) {
    return 0;
  }
  else {
    double nlt = 0;
    int n3;
    // nl term from other T
    for (int n2 = 1; n2 < NMAX; n2++) {
      n3 = n - n2;
      if (n3 > 0) {
	nlt += - PI/(2*a) * (-n2*(this->val[n3][k+1] - this->val[n3][k-1])/(2*dz) * w.val[n2][k]
			     +n3*this->val[n3][k] * (w.val[n2][k+1] - w.val[n2][k-1])/(2*dz));
      }
      n3 = n2 -n;
      if (n3 > 0) {
	nlt += -PI/(2*a)*(n2*(this->val[n3][k+1] - this->val[n3][k-1])/(2*dz) * w.val[n2][k]
			  +n3*this->val[n3][k] * (w.val[n2][k+1] - w.val[n2][k-1])/(2*dz));
      }
      n3 = n+n2;
      if (n3 > 0) {
	nlt += PI/(2*a)*(n2*(this->val[n3][k+1] - this->val[n3][k-1])/(2*dz) * w.val[n2][k]
			  +n3*this->val[n3][k] * (w.val[n2][k+1] - w.val[n2][k-1])/(2*dz));
      }
    }
    return nlt;
  }
}
  

// This function uses the default constructor of Vector
// and calculate some matrix elements
Stream::Stream (double (*init)(int, int)) :
  Vector (init)  {

  // Calculate dia[n], sup[n], sub[n]
  for (int n = 1; n < NMAX; n++){
    // Calculate each component of dia[n]
     
    double dia_val = pow(this->n * PI / a, 2) + 2*oodz_2;
    for (int k = 1; k < NZ-1; k++){
      sub[n][k] = -oodz_2;
      sup[n][k] = -oodz_2;
      dia[n][k] = dia_val;
    }
    dia[n][0   ] = 1;
    dia[n][NZ-1] = 1;
    sup[n][0   ] = 0;
    sub[n][NZ-1] = 0;
  }
    
}

void Stream::linear_step(Vector w) {
  // solve using the tridiagonal method for each mode
  for (int n = 1; n < NMAX; n++) {
    triDiSolve(w.val[n], this->val[n], this->sub[n], this->dia[n], this->sup[n]);
  }
  this->boundaries(null_0);
}

