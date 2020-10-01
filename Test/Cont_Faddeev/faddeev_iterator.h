#ifndef FADDEEV_ITERATOR_H
#define FADDEEV_ITERATOR_H

void iterate_faddeev(double** state_3N_symm_array,
                     int& Np, double** p, double** wp,
                     int& Nq, double** q, double** wq,
                     int& Nalpha, int** L_2N, int** S_2N, int** J_2N, int** T_2N, int** l_3N, int** two_j_3N);

#endif // FADDEEV_ITERATOR_H