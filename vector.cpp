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
Vector::Vector(double (*init)(int, int), double top, double bot, Simulation* dady){
  for (int n = 1; n < NMAX; n++){
    for (int k = 0; k < NZ ; k++){
      val[n][k] = init(n, k);
    }
  }
  this->top_boundary = top;
  this->bot_boundary = bot;
  this->p = dady;
}

// set boundary conditions at k = NZ-1 and k = 0
void Vector::boundaries(double (*val)(int)){
  // set boundaries
  for (int n = 1; n < NMAX; n++){
    this->val[n][NZ-1] = top_boundary;
    this->val[n][0   ] = bot_boundary;
  }

  // set mode 0
  for (int k = 0; k < NZ ; k ++)
    this->val[0][k] = val(k);
}

 
void Temp::compute_G(int n, int k){
  G_old[n][k] = G[n][k];

  G[n][k] = (this->val[n][k+1]
	     - 2*this->val[n][k]+this->val[n][k-1]) * OODZ_2 
    - pow(n * PI / a, 2) * this->val[n][k];
}

// Iterates over all modes and compute the new step
void Temp::step(){
  double nlt;
  for (int n = 0; n < NMAX; n++) {
    for (int k = 1; k < NZ-1; k++){
      // compute G
      compute_G(n,k);
	
      // Non linear term
      nlt = non_linear_term(n, k);

      // New value
      this->val[0][k] += nlt + p->dt/2 * (3*G[n][k] - G_old[n][k]);
    }
  }
  // set the boundary conditions and apply T_0 for 0th mode
  boundaries(T_0);
}

double Temp::non_linear_term(int n, int k){
  Stream* psi = p->psi;
  // nl term from T0
  double nlt = 0;
  if (n == 0) {
    for (int n2 = 1; n2 < NMAX; n2++) {
      nlt += -PI/(2*a)*n2*( ( psi->val[n2][k+1] - psi->val[n2][k-1])/(2*DZ)*this->val[n2][k]
			    + psi->val[n2][k]*(this->val[n2][k+1] - this->val[n2][k-1]/(2*DZ)));
    }
    return nlt;
  }
  else {
    int n3;
    // Contribution for n' = 0, n'' = n
    nlt = -n*PI/a * psi->val[n][k]*(this->val[0][k+1] - this->val[0][k-1])/(2*DZ);

    // nl term from other n', n''
    for (int n2 = 1; n2 < NMAX; n2++) {
      n3 = n - n2;
      if (n3 > 0 && n3 < NMAX) {
	nlt += -n*PI/a * psi->val[n][k]*(this->val[0][k+1] - this->val[0][k-1])
	  - PI/(2*a) * (-n2*(psi->val[n3][k+1] - psi->val[n3][k-1])/(2*DZ) * this->val[n2][k]
			+n3*psi->val[n3][k] * (this->val[n2][k+1] - this->val[n2][k-1])/(2*DZ));
      }
      n3 = n2 -n;
      if (n3 > 0 && n3 < NMAX) {
	nlt += -PI/(2*a)*(n2*(psi->val[n2][k+1] - psi->val[n3][k-1])/(2*DZ) * this->val[n2][k]
			  +n3*psi->val[n3][k] * (this->val[n2][k+1] - this->val[n2][k-1])/(2*DZ));
      }
      n3 = n+n2;
      if (n3 > 0 && n3 < NMAX) {
	nlt += -PI/(2*a)*(n2*(psi->val[n2][k+1] - psi->val[n3][k-1])/(2*DZ) * this->val[n2][k]
			  +n3*psi->val[n3][k] * (this->val[n2][k+1] - this->val[n2][k-1])/(2*DZ));
      }
    }
    return nlt;
  }
}

void Vort::compute_G(int n, int k){

  G_old[n] = G[n];
  G[n][k] = p->Ra*p->Pr*(n * PI / NX) * this->p->T->val[n][k] 
    + p->Pr*((this->val[n][k+1] - 2*this->val[n][k]+this->val[n][k-1]) * OODZ_2
	  - pow(n * PI / NX, 2) * this->val[n][k]);
}

  
void Vort::step(){
  double nlt;
  for (int n = 0; n < NMAX; n++) {
    for (int k = 1; k < NZ-1; k++){
      // compute G
      compute_G(n, k);
	
      // Non linear term
      nlt = non_linear_term(n, k);

      // New value
      this->val[0][k] += nlt + p->dt/2 * (3*G[n][k] - G_old[n][k]);
    }
  }
  // set the boundary conditions and apply 0 for 0th mode
  boundaries(null_0);
}

double Vort::non_linear_term(int n, int k){
  Vort* w = this->p->w;
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
	nlt += - PI/(2*a) * (-n2*(this->val[n3][k+1] - this->val[n3][k-1])/(2*DZ) * w->val[n2][k]
			     +n3*this->val[n3][k] * (w->val[n2][k+1] - w->val[n2][k-1])/(2*DZ));
      }
      n3 = n2 -n;
      if (n3 > 0) {
	nlt += -PI/(2*a)*(n2*(this->val[n3][k+1] - this->val[n3][k-1])/(2*DZ) * w->val[n2][k]
			  +n3*this->val[n3][k] * (w->val[n2][k+1] - w->val[n2][k-1])/(2*DZ));
      }
      n3 = n+n2;
      if (n3 > 0) {
	nlt += PI/(2*a)*(n2*(this->val[n3][k+1] - this->val[n3][k-1])/(2*DZ) * w->val[n2][k]
			  +n3*this->val[n3][k] * (w->val[n2][k+1] - w->val[n2][k-1])/(2*DZ));
      }
    }
    return nlt;
  }
}
  

// This function uses the default constructor of Vector
// and calculate some matrix elements
Stream::Stream (double (*init)(int, int), double top, double bot,
		Simulation* dady) :
  Vector (init, top, bot, dady)  {

  // Calculate dia[n], sup[n], sub[n]
  for (int n = 1; n < NMAX; n++){
    // Calculate each component of dia[n]
     
    double dia_val = pow(n * PI / a, 2) + 2*OODZ_2;
    for (int k = 1; k < NZ-1; k++){
      sub[n][k] = -OODZ_2;
      sup[n][k] = -OODZ_2;
      dia[n][k] = dia_val;
    }
    dia[n][0   ] = 1;
    dia[n][NZ-1] = 1;
    sup[n][0   ] = 0;
    sub[n][NZ-1] = 0;
  }
    
}

void Stream::step() {
  // solve using the tridiagonal method for each mode
  for (int n = 1; n < NMAX; n++) {
    triDiSolve(this->p->w->val[n], this->val[n], this->sub[n], this->dia[n], this->sup[n]);
  }
  boundaries(null_0);
}
