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
  
  void step();
  void boundaries(double (*val)(int));

  Vector(double (*f)(int, int), double top, double bot,
	 Simulation* dady);
};

// Describes a temperature field
class Temp : public Vector{
  using Vector::Vector;
private:
  double non_linear_term(int, int);
  void compute_G(int, int);
public:
  void step();

};

// Describes a vorticity field
class Vort : public Vector{
  using Vector::Vector;
private:
  double non_linear_term(int, int);
  void compute_G(int, int);
public:
  void step();
};

// Describes a stream field
class Stream : public Vector{
  using Vector::Vector;
private:
  field sub, sup, dia;
public:
  void step();

  Stream(double (*init)(int, int), double, double, Simulation*);
};

class Simulation {
private:
  std::array<double,NMAX> old_T, old_w, old_psi;
  void dump();
  void cfl();

public:
  Temp* T;
  Vort* w;
  Stream* psi;
  int nmax, niter, save_freq, save_freq_growth, max_niter;
  double max_time, dt, t, dt_security;
  // Constants
  double Pr, Ra, Re;
  // Iteration method, returns true as long as we want
  // and save the required stuff
  bool iter();

  Simulation(int save_freq = 500, int save_freq_growth = 0,
	     int max_niter = 10000, double max_time = 0,
             double init_time = 0, int init_niter = 0,
             double dt_security = 0.9);
  ~Simulation() {
    delete T;
    delete w;
    delete psi;
  }
};
    
#endif
