#include "misc.hpp"

// Describes the nth component of the x-FFT of a 2-D field
class Vector{
protected:
  field val;
  double top_boundary, bot_boundary;
  
public:
  field G, G_old;
  int n;

  Vector(double (*f)(int, int), double top = 0, double bot = 0);
  vector mode(int n) {return val[n];};
  void boundaries(double (*val)(int));
};

// Describes a temperature field
class Temp : public Vector{
  using Vector::Vector;
private:
  void compute_G(int);
public:
  void step();
};

// Describes a vorticity field
class Vort : public Vector{
  using Vector::Vector;
private:
  void compute_G(int,Temp);
public:
  void step(Temp);
};

// Describes a stream field
class Stream : public Vector{
  using Vector::Vector;
public:
  field sub, sup, dia;

  void step(Vort);
  Stream(double (*init)(int, int));
};

    
class Simulation{
private:
  std::array<double,NMAX> old_T, old_w, old_phi;
  void dump(Temp, Vort, Stream );
public:
  double time = 0;
  int nmax; // Max precision in fourier-space
  int niter;
  int save_freq;
  // Iteration method, returns true as long as we want
  // and save the required stuff
  bool iter(Temp, Vort, Stream);

  Simulation(int save_freq = 10, double init_time = 0,
			int init_niter = 0);
};
