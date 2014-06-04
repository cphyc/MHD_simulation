#ifndef MISC_H
#define MISC_H

#include <array>
#include "constant.cpp"

// Type of our vectors (!= std::vector)
typedef std::array<double,NZ> vector;
// Type of our fields (array of vectors)
typedef std::array<vector,NMAX> field;

void triDiSolve(vector& rhs, vector &sol, vector& sub,
		vector& dia, vector& sup);

// void null_vector();

// Initialization functions

double T_init (int n, int k);
double T_0 (int k);

double null_init(int n, int k);
double null_0(int k);

double abs(double x) ;

void pp(double) ;
#endif
