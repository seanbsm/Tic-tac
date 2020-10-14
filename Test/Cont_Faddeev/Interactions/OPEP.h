#ifndef OPEP_H
#define OPEP_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "../constants.h"
#include "../General_functions/legendre.h"
#include "../General_functions/gauss_legendre.h"

class OPEP
{
private:
	int nz = 150;	// number of points used in numerical angular integral
	
	std::vector<double> z;
	std::vector<double> wz;

	LegendrePolynomial* legPol = new LegendrePolynomial();

	int    isoFac(int L, int S);
	double angIntegral(double qi, double qo, int J, int l);
	double potOPEPmom(double qi, double qo, double z);
	
public:
	OPEP(int Jmin, int Jmax);
	
	void potential(double qi, double qo, int J, double *Varray);
};

#endif // OPEP_H

