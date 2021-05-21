#ifndef KINETIC_CONVERSION_H
#define KINETIC_CONVERSION_H

#include "../type_defs.h"
#include "../constants.h"

double com_momentum_to_lab_energy(double q_com);
double lab_energy_to_com_momentum(double T_lab);

double lab_energy_to_com_energy(double T_lab);
double com_energy_to_lab_energy(double T_com);

double com_energy_to_com_q_momentum(double E_com);
double com_q_momentum_to_com_energy(double q_com);

#endif // KINETIC_CONVERSION_H