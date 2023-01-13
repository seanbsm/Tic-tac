
#include "kinetic_conversion.h"

/* Deuteron mass (depends on binding energy) */
double deuteron_mass(double Eb){
	return Mn + Mp + Eb; // Deuteron mass [MeV]
}

/* Reduced mass for deuteron-nucleon c.m. system */
double mu1(double Eb){
	double Md = deuteron_mass(Eb);
	return Mn*Md/(Mn + Md);
}

/* Reduced mass for nucleon-nucleon-nucleon c.m. system */
double mu1_free(){
	double M12 = Mp+Mn;
	return Mn*M12/(Mn + M12);
}

/* Converts CoM momentum to Lab kinetic energy */
double com_momentum_to_lab_energy(double q_com, double Eb){
	
	double Md = deuteron_mass(Eb);

	double b = 2*(Mp - q_com*q_com/Md);
	double c = -q_com*q_com*(Mp+Md)*(Mp+Md)/(Md*Md);
	
	return 0.5*(-b + sqrt(b*b - 4*c));
}

/* Converts Lab kinetic energy to CoM momentum (assuming p+d scattering) */
double lab_energy_to_com_momentum(double T_lab, double Eb){
	double Md = deuteron_mass(Eb);
	return sqrt( Md*Md*T_lab*(T_lab + 2*Mp) / ( (Mp+Md)*(Mp+Md) + 2*T_lab*Md ) );
}

/* Converts 3N CoM kinetic energy to q-momentum for nucleon-deuteron scattering */
double com_energy_to_com_q_momentum(double E_com){
	return sqrt(4*E_com*MN/3);
}

/* Converts 3N q-momentum to CoM kinetic energy for nucleon-deuteron scattering */
double com_q_momentum_to_com_energy(double q_com){
	return 3*q_com*q_com/(4*MN);
}


double com_energy_to_com_p_momentum(double Ep){
	return sqrt(2*Ep*mu23);
}
double com_p_momentum_to_com_energy(double p){
	return p*p/(2*mu23);
}

//double lab_energy_to_com_energy(double T_lab, int Tz){
//}
//
//double com_energy_to_lab_energy(double T_com, int Tz){
//}
