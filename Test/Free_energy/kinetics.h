#ifndef KINETICS_H
#define KINETICS_H

void get_p_and_q_momenta (double* p_array,
                          double* q_array,
                          int n);

void calculate_3N_kinetic_energy(double* p_array,
                                 double* q_array,
                                 double* T_array,
                                 int n);

#endif // KINETICS_H