
#include "kinetics.h"

void get_p_and_q_momenta (double* p_array,
                          double* q_array,
                          int n){

}

void calculate_3N_kinetic_energy(double* p_array,
                                 double* q_array,
                                 double* T_array,
                                 int n){
    
    double m1 = 1;
    double m2 = 1;
    double m3 = 1;
    
    double mu1 = 1/m1 + 1/m2;
    double mu2 = 1/(m1+m2) + 1/m3;

}