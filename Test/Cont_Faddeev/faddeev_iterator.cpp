
#include "faddeev_iterator.h"

double p1(double q, double qp, double x){
    return sqrt(0.25*q*q + qp*qp + x*q*qp); 
}

double p2(double q, double qp, double x){
    return sqrt(q*q + 0.25*qp*qp + x*q*qp);
}

void iterate_faddeev(double** state_3N_symm_array,
                     int& Np, double** p, double** wp,
                     int& Nq, double** q, double** wq,
                     int& Nalpha, int** L_2N, int** S_2N, int** J_2N, int** T_2N, int** l_3N, int** two_j_3N){
    
    double E;
    double pre_fac = 1./ (E - p*p/(2*mu) - q*q/(2*M));
}