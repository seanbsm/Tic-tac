#ifndef READ_PSI_SEAN_H
#define READ_PSI_SEAN_H

#include <iostream>
#include <stdio.h>

#include "../error_management.h"
#include "../General_functions/spin_coupling_functions.h"

void read_psi (char* filename, 
               double** psi_read,
               int *Np, double **p, double **wp,
               int *Nq, double **q, double **wq,
               int *Nalpha, int **get_L, int **get_S, int **get_J, int **get_T, int **get_l, int **get_two_j);

void get_all_states(double** state_3N_symm_array,
                    double** state_3N_asym_array,
                    int& Np, double** p, double** wp,
                    int& Nq, double** q, double** wq,
                    int& Nalpha, int** L_2N, int** S_2N, int** J_2N, int** T_2N, int** l_3N, int** two_j_3N);

#endif // READ_PSI_SEAN_H