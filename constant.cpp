#ifndef CONST
#define CONST

#include <cmath>

inline double max(double a, double b){ if (a > b) return a; else return b;}

const int   NZ        =         100; // Z dimension
const int   a         =           1;
const int   NX        =        a*NZ; // X dimension
const int   NMAX      =          50; // Resolution for FFT

const double DZ       =  1.0/(NZ-1); // Vertical step
const double OODZ_2   = 1.0/(DZ*DZ);// 1/dz²
const double PI       = 4.0*atan(1);
const double c        =        PI/a;
#endif
