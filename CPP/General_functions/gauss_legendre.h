#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include <iostream>
#include <vector>
#include <cmath>

void gauss( float* x,
			float* w,
			int N );
void gauss( double* x,
			double* w,
			int N );

void rangeChange_0_inf( float* x,
						float* w,
						float scale,
						int N );
void rangeChange_0_inf( double* x,
						double* w,
						double scale,
						int N );

void updateRange_a_b( float* x,
					  float* w,
					  float a,
					  float b,
					  int N );
void updateRange_a_b( double* x,
					  double* w,
					  double a,
					  double b,
					  int N );

#endif // GAUSS_LEGENDRE_H

