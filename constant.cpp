#ifndef CONST
#define CONST

#include <cmath>

inline double max(double a, double b){ if (a > b) return a; else return b;}

const int   NZ        =         100; // Z dimension
const int   a         =           1;
const int   NX        =        a*NZ; // X dimension
const int   NMAX      =          50; // Resolution for FFT

extern double Ra, Pr; // declare Ra and Pr

const double dz       =  1.0/(NZ-1); // Vertical step
const double oodz_2   = 1.0/(dz*dz);// 1/dz²
extern double dt;
const double MAX_TIME =        1000; // Max time
extern int   MAX_NITER;
const double PI       = 4.0*atan(1);
const int FREQ_GROWTH =         500;
const double c        =        PI/a;

const double DT       =         1.0;
#endif
