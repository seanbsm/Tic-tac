#ifndef CHIRAL_LO_SIMPLE_H
#define CHIRAL_LO_SIMPLE_H

#include <iostream>
#include <vector>
#include <cmath>

class chiral_LO
{
private:

    /* Gauss-Legendre grid objects for angular OPEP integral */
    int nz  = 150;
    std::vector<double> z;
	std::vector<double> wz;
    /* Vectors holding the roots of the Legendre polynomials */
    std::vector<std::vector<double>> roots_d;

    /* Natural constants */
    const double Mp		= 938.272;			// proton mass  						[MeV]
    const double Mn		= 939.565;			// neutron mass 						[MeV]
    const double MN		= 2*Mp*Mn/(Mp+Mn);	// Nucleon mass							[MeV]
    const double gA 	= 1.289;			// axial coupling constant 			 	[no units]
    const double fpi 	= 92.2;				// pion decay constant 				 	[MeV] 		(used convention 130.41/sqrt(2) = 92.2)
    const double mpi 	= 138.039;			// pion_0 mass 						 	[MeV] 		(average of pi+, pi-, and pi0 masses)
    const double hbarc  = 197.326;          // MeV fm
    /* Chiral LO parameters */
    const double Lambda	= 450;				// cut-off for renormalization of LO 	[MeV]
    /* Chiral LO LECs */
    const double C1S0	= -0.112927/100.;	// contact term C1S0 for lambda = 450	[MeV]
    const double C3S1	= -0.087340/100.;	// contact term C3S1 for lambda = 450	[MeV]

    /* Functions to construct vectors with roots of Legendre polynomials */
    void find_roots(int L, std::vector<double> &z);
    void fetch_root(int d, int i, double &root);
    /* Function to construct a Gauss-Legendre grid from x=-1 to x=1, with weigths w, of size N */
    void make_gauss_legendre_grid(double* x, double* w, int N);

    double calculate_nucleon_mass(int& Tz);

    /* Functions for calculating one-pion exchange potential */
    int isospin_prefactor(int L, int S);
    void OPEP_potential(double qi, double qo, int J, double *Varray);
    double OPEP_angular_integral(double qi, double qo, int J, int l);
    double OPEP_momentum_terms(double qi, double qo, double z);
public:
	chiral_LO();
	
	void V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int& Tz, double *Varray);
};

#endif // CHIRAL_LO_SIMPLE_H
