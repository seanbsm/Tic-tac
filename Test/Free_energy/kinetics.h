#ifndef KINETICS_H
#define KINETICS_H

const double hbarc = 197.326; // MeV fm
const double MN = 938.918; // in MeV, averaged mass (Mn + Mp)/2

double calculate_3N_kinetic_energy(double* state_3N_symm_array,
                                   double* state_3N_asym_array,
                                   int &Np, double* p_array, double* wp_array,
                                   int &Nq, double* q_array, double* wq_array,
                                   int& Nalpha);

#endif // KINETICS_H