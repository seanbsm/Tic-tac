#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <iostream>
#include <vector>
#include <map>

#include "../typedefs.h"

class LegendrePolynomial
{
private:
	 std::vector<std::vector<float>>  roots_f;
	 std::vector<std::vector<double>> roots_d;
public:
	LegendrePolynomial();
	
	/* Calculates a single root */
	void findRoot(int d, float x, float &root);
	void findRoot(int d, double x, double &root);
	
	/* Calculates a range of roots, one for each z[i] and
	 * l in [L1,L2] */
	void findRoots(int L, std::vector<float> &z);
	void findRoots(int L, std::vector<double> &z);
	
	/* Retrieves a single root, precalculated by findRoots */
	void fetchRoot(int d, int i, float &root);
	void fetchRoot(int d, int i, double &root);
};
#endif // LEGENDRE_H

