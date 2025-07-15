
#include "OPEP.h"

OPEP::OPEP(int Jmin, int Jmax){
	
	z.resize(nz);
	wz.resize(nz);
	
	/* Find a Gauss-Legendre distribution */
	gauss(&z[0],&wz[0],nz);

	/* Precalculate all legendre polynomials
	 * required for angular integral */
	legPol->findRoots(Jmax+1, z);
}

int OPEP::isoFac(int L, int S){
	int T = ( (1-L-S) & 1);
	int isoFactor = -3*(1-T) + 1*T;
	
	return isoFactor;
}

void OPEP::potential(double qi, double qo, int J, double *Varray){
	
	double integral_0 = angIntegral(qi,qo,J,0);
	double integral_1 = angIntegral(qi,qo,J,1);
	double integral_P = angIntegral(qi,qo,J+1,0);
	
	double V_uncoupled_S0 = 0;
	double V_uncoupled_S1 = 0;
	double V_coupled_mm   = 0;
	double V_coupled_pm   = 0;
	double V_coupled_mp   = 0;
	double V_coupled_pp   = 0;
	
	V_uncoupled_S0 = 2 * (-(qo*qo+qi*qi)*integral_0 + 2*qo*qi*integral_1);
	V_coupled_pp   = (2./(2*J+1)) * (-(qo*qo+qi*qi)*integral_P + 2*qo*qi*integral_0);
	
	V_uncoupled_S0 *= isoFac(J,0);
	V_coupled_pp   *= isoFac(J+1,1);
	
	if (J!=0){
		double integral_M = angIntegral(qi,qo,J-1,0);
		
		V_uncoupled_S1 = 2 * ((qo*qo+qi*qi)*integral_0 - 2*qo*qi*(1./(2*J+1))*(J*integral_P + (J+1)*integral_M));
		V_coupled_mm = (2./(2*J+1)) * ((qo*qo+qi*qi)*integral_M - 2*qo*qi*integral_0);
		V_coupled_pm = (4*sqrt(J*(J+1))/(2*J+1)) * (qi*qi*integral_P + qo*qo*integral_M - 2*qo*qi*integral_0);
		V_coupled_mp = (4*sqrt(J*(J+1))/(2*J+1)) * (qi*qi*integral_M + qo*qo*integral_P - 2*qo*qi*integral_0);
		
		V_uncoupled_S1 *= isoFac(J,1);
		V_coupled_mm *= isoFac(J-1,1);
		V_coupled_pm *= isoFac(J-1,1);
		V_coupled_mp *= isoFac(J+1,1);
	}
	
	Varray[0] = V_uncoupled_S0;
	Varray[1] = V_uncoupled_S1;
	Varray[2] = V_coupled_mm;
	Varray[3] = V_coupled_pm;
	Varray[4] = V_coupled_mp;
	Varray[5] = V_coupled_pp;
}
	
double OPEP::angIntegral(double qi, double qo, int J, int l){
	
	double integral = 0;
	double root = 0;
	if (l==0){
		for (size_t i=0; i<wz.size(); i++){
			legPol->fetchRoot(J,i,root);
			integral += wz[i] * potOPEPmom(qi,qo,z[i]) * root;
		}
	}
	else{
		for (size_t i=0; i<wz.size(); i++){
			legPol->fetchRoot(J,i,root);
			integral += wz[i] * potOPEPmom(qi,qo,z[i]) * z[i] * root;
		}
	}
	
	return pi*integral;
}

double OPEP::potOPEPmom(double qi, double qo, double z){
	double q2 = qi*qi + qo*qo - 2*qi*qo*z;
	
	return -(gA*gA/(4*fpi*fpi))*(1./(q2 + mpi*mpi));
}


