#ifndef vector_h
#define vector_h

#include <fftw3.h>
#include "misc.hpp"

class Simulation;

// Describes the nth component of the x-FFT of a 2-D field
class Vector{
protected:
  Simulation* p;
  double top_boundary, bot_boundary;
  field nonlinear_contribution();
  double in[NX], out[NX];
  fftw_r2r_kind fft_flag;
  bool FFT_has_been_done;
  
public:

  void unset_FFT(){ FFT_has_been_done = false; }
  field val, dval;
  field G, G_old;
  std::array<std::array<double,NZ>, NX> real_val;
  int n;
  
  void compute();
  void add();
  void boundaries();
  void FFT();

  Vector(double top, double bot,
	 Simulation* dady);
};

// Describes a temperature field
class Temp : public Vector{
private:
  double non_linear_term(int, int);
  void compute_G(int, int);
public:
  void compute();
  Temp (double top, double bot, Simulation* dady) : 
    Vector(top, bot, dady)
  { this->fft_flag = FFTW_REDFT00; } 
};

// Describes a vorticity field
class Vort : public Vector{
private:
  double non_linear_term(int, int);
  void compute_G(int, int);
public:
  void compute();
  Vort (double top, double bot, Simulation* dady) : 
    Vector(top, bot, dady)
  { this->fft_flag = FFTW_RODFT00; } 
};

// Describes a stream field
class Stream : public Vector{
  using Vector::Vector;
private:
  int fft_flag = FFTW_RODFT00;
  field sub, sup, dia;
public:
  void compute();

  Stream(double, double, Simulation*);
};

class Simulation {
private:
  std::array<double,NMAX> old_T, old_w, old_psi;
  void dump();
  void cfl();
  void unset_FFT();
  void do_FFT();

public:
  Temp* T;
  Vort* w;
  Stream* psi;
  int nmax, niter, save_freq, save_freq_growth, max_niter, cfl_freq;
  double max_time, dt, dt_old, t, dt_security, dt_additional_cfl_security;
  // Constants
  double Pr, Ra, Re;
  // Iteration method, returns true as long as we want
  // and save the required stuff
  bool iter();

  Simulation(int save_freq = 500, int save_freq_growth = 0,
	     int max_niter = 10000, double max_time = 0,
             double init_time = 0, int init_niter = 0,
             double dt_security = 0.8, double dt_additional_cfl_security = 0.25,
	     double init_dt = 1e-6);
  ~Simulation() {
    delete T;
    delete w;
    delete psi;
  }
};
    
#endif
