
#include "chiral_LO_simple.h"

chiral_LO::chiral_LO(){
	/* Make a Gauss-Legendre distribution */
	make_gauss_legendre_grid(&z[0],&wz[0],nz);
    z.resize(nz);
	wz.resize(nz);

    /* Precalculate all legendre polynomials
	 * required for angular integral */
	find_roots(100, z);
}

void chiral_LO::find_roots(int L, std::vector<double> &z){
	
	/* List with roots per l, for each z[i] */
	std::vector<double> currRootList; (z.size());
	
	/* Set all roots to 1 for l=0 */
	currRootList.assign(z.size(), 1);
	roots_d.push_back(currRootList);
	
	/* Set all roots to z[i] for l=1 */
	roots_d.push_back(z);
	
	/* Use recurrance relation for Legendre polynomials
	 * to set remainder of roots */
	for (int l=1; l<L; l++){
		for (size_t i=0; i<z.size(); i++){
			double Pcurr = roots_d[l][i];
			double Pprev = roots_d[l-1][i];
			double root = ( (2*l+1)*z[i]*Pcurr - l*Pprev ) / (l+1);
			currRootList[i] = root;
		}
		roots_d.push_back(currRootList);
	}
}

void chiral_LO::fetch_root(int d, int i, double &root){
	root = roots_d[d][i];
}

void chiral_LO::make_gauss_legendre_grid(double* x, double* w, int N){
	double p3	= 0;								// Legendre polynomial, P(xi), for i=3
	double eps = 3.e-16;
	int m	= (N+1)/2;
	
	for (int i=1; i<m+1; i++){
		double t  = cos(M_PI*(i-0.25)/(N+0.5));		// xi
		double t1 = 1;								// old xi
		double pp = 0;
		
		while (std::abs(t-t1) >= eps){
			double p1 = 1;							// Legendre polynomial, P(xi), for i=1
			double p2 = 0;							// Legendre polynomial, P(xi), for i=2
			
			for (int j=1; j<N+1; j++){
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*t*p2 - (j-1)*p3)/j;	// recurrence relation for Legendre polynomials
			}
			pp = N*(t*p1-p2)/(t*t-1);				// identity for P'(xi)
			t1 = t;
			t  = t1 - p1/pp; 						// Newton's method for finding roots
		}
		
		x[i-1] = -t;
		x[N-i] = t;
		
		w[i-1] = 2./((1-t*t)*pp*pp);
		w[N-i] = w[i-1];
	}
}

double chiral_LO::calculate_nucleon_mass(int& Tz){
    /* Nucleon mass */
    double MN = 0;

    /* Define Nucleon mass (as defined by Machleidt) */
    if (Tz==0){
        MN = 2*Mn*Mp/(Mn+Mp);
    }
    else if (Tz==-1){
        MN = Mp;
    }
    else if (Tz==+1){
        MN = Mn;
    }
    else{
        std::cout << "Encountered invalid value Tz=" << Tz << " in chiral_LO nucleon-mass calculation." << std::endl;
        std::cout << "Returning nucleon mass MN=0 MeV" << std::endl;
    }

	return MN;
}

void chiral_LO::V(double& qi, double& qo, bool coupled, int& S, int& J, int& T, int& Tz, double* Varray){
    
	Varray[0] = 0;
	Varray[1] = 0;
	Varray[2] = 0;
	Varray[3] = 0;
	Varray[4] = 0;
	Varray[5] = 0;
	
	/* Pion-exhange ("OPEP" = One-Pion Exchange Potential) */
	OPEP_potential(qi, qo, J, Varray);
	
	/* Contact terms */
	if (J==0){
		Varray[0] += C1S0;
	}
	else if (J==1){
		Varray[2] += C3S1;
	}

    /* Nucleon mass */
    double MN = calculate_nucleon_mass(Tz);
	
    /* Coefficients variable */
    double coeff = 1;

	/* Minimal relativity factors */
	double Epi         = sqrt(MN*MN + qi*qi); 	// relativistic energy of in-going particle
	double Epo         = sqrt(MN*MN + qo*qo); 	// relativistic energy of out-going particle
	double relFactor_i = sqrt(MN/Epi);      	// relativistic factor of in-going particle
	double relFactor_o = sqrt(MN/Epo);      	// relativistic factor of out-going particle
	coeff *= relFactor_i*relFactor_o;

	/* Regulator-functions */
	double temp1 = qi/Lambda;
	double temp2 = temp1*temp1*temp1*temp1*temp1*temp1;	// i.e. we have (qi/Lambda)^6
	double f1 	 = exp(-temp2);

	temp1        = qo/Lambda;
	temp2        = temp1*temp1*temp1*temp1*temp1*temp1;	// i.e. we have (qo/Lambda)^6
	double f2    = exp(-temp2);
	coeff *= f1*f2;
    
    /* Fourier transform constants */
    coeff /= (8*M_PI*M_PI*M_PI);

	Varray[0] *= coeff;
	Varray[1] *= coeff;
	Varray[2] *= coeff;
	Varray[3] *= coeff;
	Varray[4] *= coeff;
	Varray[5] *= coeff;
}

int chiral_LO::isospin_prefactor(int L, int S){
	int T = ( (1-L-S) & 1);
	int isoFactor = -3*(1-T) + 1*T;
	
	return isoFactor;
}

void chiral_LO::OPEP_potential(double qi, double qo, int J, double *Varray){
	
	double integral_0 = OPEP_angular_integral(qi,qo,J,0);
	double integral_1 = OPEP_angular_integral(qi,qo,J,1);
	double integral_P = OPEP_angular_integral(qi,qo,J+1,0);
	
	double V_uncoupled_S0 = 0;
	double V_uncoupled_S1 = 0;
	double V_coupled_mm   = 0;
	double V_coupled_pm   = 0;
	double V_coupled_mp   = 0;
	double V_coupled_pp   = 0;
	
	V_uncoupled_S0 = 2 * (-(qo*qo+qi*qi)*integral_0 + 2*qo*qi*integral_1);
	V_coupled_pp   = (2./(2*J+1)) * (-(qo*qo+qi*qi)*integral_P + 2*qo*qi*integral_0);
	
	V_uncoupled_S0 *= isospin_prefactor(J,   0);
	V_coupled_pp   *= isospin_prefactor(J+1, 1);
	
	if (J!=0){
		double integral_M = OPEP_angular_integral(qi,qo,J-1,0);
		
		V_uncoupled_S1 = 2 * ((qo*qo+qi*qi)*integral_0 - 2*qo*qi*(1./(2*J+1))*(J*integral_P + (J+1)*integral_M));
		V_coupled_mm   = (2./(2*J+1)) * ((qo*qo+qi*qi)*integral_M - 2*qo*qi*integral_0);
		V_coupled_pm   = (4*sqrt(J*(J+1))/(2*J+1)) * (qi*qi*integral_P + qo*qo*integral_M - 2*qo*qi*integral_0);
		V_coupled_mp   = (4*sqrt(J*(J+1))/(2*J+1)) * (qi*qi*integral_M + qo*qo*integral_P - 2*qo*qi*integral_0);
		
		V_uncoupled_S1 *= isospin_prefactor(J,   1);
		V_coupled_mm   *= isospin_prefactor(J-1, 1);
		V_coupled_pm   *= isospin_prefactor(J-1, 1);
		V_coupled_mp   *= isospin_prefactor(J+1, 1);
	}
	
	Varray[0] = V_uncoupled_S0;
	Varray[1] = V_uncoupled_S1;
	Varray[2] = V_coupled_mm;
	Varray[3] = V_coupled_pm;
	Varray[4] = V_coupled_mp;
	Varray[5] = V_coupled_pp;
}
	
double chiral_LO::OPEP_angular_integral(double qi, double qo, int J, int l){
	
	double integral = 0;
	double root = 0;
	if (l==0){
		for (size_t i=0; i<wz.size(); i++){
			fetch_root(J,i,root);
			integral += wz[i] * OPEP_momentum_terms(qi,qo,z[i]) * root;
		}
	}
	else{
		for (size_t i=0; i<wz.size(); i++){
			fetch_root(J,i,root);
			integral += wz[i] * OPEP_momentum_terms(qi,qo,z[i]) * z[i] * root;
		}
	}
	
	return M_PI*integral;
}

double chiral_LO::OPEP_momentum_terms(double qi, double qo, double z){
	double q2 = qi*qi + qo*qo - 2*qi*qo*z;
	
	return -(gA*gA/(4*fpi*fpi))*(1./(q2 + mpi*mpi));
}


