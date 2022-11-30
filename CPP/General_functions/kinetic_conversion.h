#ifndef KINETIC_CONVERSION_H
#define KINETIC_CONVERSION_H

#include "../type_defs.h"
#include "../constants.h"

double deuteron_mass(double Eb);
double mu1(double Eb);

double com_momentum_to_lab_energy(double q_com, double Eb);
double lab_energy_to_com_momentum(double T_lab, double Eb);

double lab_energy_to_com_energy(double T_lab, double Eb);
double com_energy_to_lab_energy(double T_com, double Eb);

double com_energy_to_com_q_momentum(double E_com);
double com_q_momentum_to_com_energy(double q_com);

double com_energy_to_com_p_momentum(double Ep);
double com_p_momentum_to_com_energy(double p);

#endif // KINETIC_CONVERSION_H