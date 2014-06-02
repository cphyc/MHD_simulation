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

 
void Temp::compute_G(int n){
     G_old[n] = G[n];
     for (int k = 1; k < NZ-1; k++){
	  G[n][k] = (val[n][k+1]
		     - 2*val[n][k]+val[n][k-1]) * oodz_2 
	       - pow(n * PI / a, 2) * val[n][k];
     }
}

// Iterates over all modes and compute the new step
void Temp::step(){
     for (int n = 1; n < NMAX; n++) {
	  // compute G_n(t) and G_n(t-dt)
	  compute_G(n);
	  // calculate the new T
	  for (int k = 0; k < NZ; k++){
	       this->val[n][k] = 
		    this->val[n][k] + dt/2 * (3*G[n][k] - G_old[n][k]);
	       //std::cout << "Out" << dz << " "<< oodz_2 << " "<< dt << std::endl;
	  }
	  // set the boundary conditions and apply T_0 for 0th mode
	  this->boundaries(T_0);
     }
}

void Vort::compute_G(int n, Temp T){
     //int k = 1;
     //std::cout << "#Duh!"<<n << " "<< T.mode(n)[k] << std::endl;
     G_old[n] = G[n];
     for (int k = 1; k < NZ-1; k++){
	  G[n][k] = Ra*Pr*(n * PI / NX) * T.val[n][k] 
	       + Pr*((val[n][k+1] - 2*val[n][k]+val[n][k-1]) * oodz_2
		     - pow(n * PI / NX, 2) * val[n][k]);
     }
}

  
void Vort::step(Temp T){
     for (int n = 1; n < NMAX; n++){
	  // compute G(t) and G(t-dt)
	  compute_G(n, T);
	  // calculate the new w
	  for (int k = 0; k < NZ; k++)
	       val[n][k] = val[n][k] + dt/2*(3*G[n][k] - G_old[n][k]);
     }
     this->boundaries(null_0);
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

void Stream::step(Vort w) {
     // solve using the tridiagonal method for each mode
     for (int n = 1; n < NMAX; n++) {
	  triDiSolve(w.val[n], this->val[n], this->sub[n], this->dia[n], this->sup[n]);
     }
     this->boundaries(null_0);
}

Simulation::Simulation (int save_freq, double init_time,
			int init_niter){
     this->save_freq = save_freq ;
     this->time = init_time; 
     this->niter = init_niter;
}


// Get an array of T, w, phi and return a string containing :
// T0 w0 phi0 T1 w1 phi1 ... Tn wn phin
void Simulation::dump(Temp T, Vort w, Stream phi) {
     // write the header
     std::cout << "\n\n#Header t:\t" << time << "\n" ;

     // count the columns
     std::cout << std::endl << "#";
     for (int n = 1; n < NMAX; n++)
	  std::cout << n << "\t";
     std::cout << std::endl;

     // print the title
     std::cout << "k" ;
     for (int n = 1; n < NMAX; n++)
	  std::cout << "\tT_" << n << "\tw_" << n << "\tphi_" << n ;
  
     // get back to the next line
     std::cout << "\n";

     // write the data (one line = one k)
     for (int k = 0; k < NZ; k++){
	  std::cout << k << "\t";
    
	  for (int n = 0; n < NMAX; n++){
	       //std::cout << "#k:n :\t" << k <<":"<< n << "\n";
	       std::cout << (T.val[n][k]) << "\t" ;
	       std::cout << (w.val[n][k]) << "\t";
	       std::cout << (phi.val[n][k]) << "\t";
	  }
	  // ended one n, we get back to the next line for the next k
	  std::cout<< "\n";
     }
}

// Iterator of the simulation. Stops when the MAX_TIME or MAX_ITER is reached.
bool Simulation::iter(Temp T, Vort w, Stream phi){  
     bool cont = true;
     // Save every save_freq
     if (niter % save_freq == 0) {
	  dump(T,w,phi);
     }
     // Save the growth rate every FREQ_GROWTH
     if (niter % FREQ_GROWTH == 0) {
	  // Growth rate calculated at a third of the heigth
	  int third = NZ / 3;

	  // Save the initial fields
	  if (niter == 0) {
	       for (int n = 1; n < NMAX; n++){
		    old_T[n] = T.val[n][third];
		    old_w[n] = w.val[n][third];
		    old_phi[n] = phi.val[n][third];
	       }
	  }
	  
	  // Iterate over all the modes
	  for (int n = 1; n < NMAX; n++){
	       // compute log | T_old | / | T |
	       double gr_T = log( abs(T.val[n][third])) - log(abs(old_T[n]) );
	       double gr_w = log(abs(w.val[n][third])) - log(abs(old_w[n]));
	       double gr_phi = log(abs(phi.val[n][third])) - log(abs(old_phi[n]));

	       // output it
	       std::cout << "# Growth rate :\t" << n << "\t";
	       std::cout << niter << "\t" << gr_T;
	       std::cout << "\t" << gr_w << "\t" << gr_phi << "\n";

	       // save it for the next output
	       old_T[n] = T.val[n][third];
	       old_w[n] = w.val[n][third];
	       old_phi[n] = phi.val[n][third];      
	  }
     }
  
     // move of one time step and check we can continue
     niter ++;
     time ++;
     if ((time > MAX_TIME) || (niter > MAX_NITER)) cont = false;
     return cont;
}
