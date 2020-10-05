#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include <iostream>
#include <vector>
#include <cmath>

void gauss( std::vector<float>& x,
			std::vector<float>& w,
			int N );
void gauss( std::vector<double>& x,
			std::vector<double>& w,
			int N );

void rangeChange_0_inf( std::vector<float>& x,
						std::vector<float>& w,
						float scale );
void rangeChange_0_inf( std::vector<double>& x,
						std::vector<double>& w,
						double scale );

void updateRange_a_b( std::vector<float>& x,
					  std::vector<float>& w,
					  float a,
					  float b );
void updateRange_a_b( std::vector<double>& x,
					  std::vector<double>& w,
					  double a,
					  double b );

#endif // GAUSS_LEGENDRE_H

