#ifndef vector_h
#define vector_h

#include "misc.hpp"

class Simulation;

// Describes the nth component of the x-FFT of a 2-D field
class Vector{
protected:
  Simulation* p;
  double top_boundary, bot_boundary;
  field nonlinear_contribution();  
public:
  field val;
  field G, G_old;
  int n;

  Vector(double (*f)(int, int), double top = 0, double bot = 0);
  void boundaries(double (*val)(int));
};

// Describes a temperature field
class Temp : public Vector{
  using Vector::Vector;
private:
  double non_linear_term(Vector, int, int);
  void compute_G(int, int);
public:
  void linear_step(Vector);

};

// Describes a vorticity field
class Vort : public Vector{
  using Vector::Vector;
private:
  double non_linear_term(Vector, int, int);
  void compute_G(Vector, int, int);
public:
  void linear_step(Vector, Vector);
};

// Describes a stream field
class Stream : public Vector{
  using Vector::Vector;
public:
  field sub, sup, dia;

  void linear_step(Vector);
  Stream(double (*init)(int, int));
};

class Simulation {
private:
  Temp *T;
  Vort *w;
  Stream *psi;
  std::array<double,NMAX> old_T, old_w, old_psi;
  void dump(Temp, Vort, Stream );
  int max_niter;
  double max_time;
public:
  double time = 0;
  int nmax; // Max precision in fourier-space
  int niter;
  int save_freq, save_freq_growth;
  // Iteration method, returns true as long as we want
  // and save the required stuff
  bool iter(Temp, Vort, Stream);

  Simulation(int save_freq = 500, int save_freq_growth = 0,
	     int max_niter = 10000, double max_time = 0,
	     double init_time = 0, int init_niter = 0);
};
    
#endif
