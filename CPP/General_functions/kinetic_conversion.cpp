
#include "kinetic_conversion.h"

/* Converts CoM momentum to Lab kinetic energy */
double com_momentum_to_lab_energy(double q_com){
	
	double b = 2*(Mp - q_com*q_com/Md);
	double c = -q_com*q_com*(Mp+Md)*(Mp+Md)/(Md*Md);
	
	return 0.5*(-b + sqrt(b*b - 4*c));
}

/* Converts Lab kinetic energy to CoM momentum (assuming p+d scattering) */
double lab_energy_to_com_momentum(double T_lab){
	return sqrt( Md*Md*T_lab*(T_lab + 2*Mp) / ( (Mp+Md)*(Mp+Md) + 2*T_lab*Md ) );
}

double lab_energy_to_com_energy(double T_lab, int Tz){
}

double com_energy_to_lab_energy(double T_com, int Tz){
}
