#include "make_wp_states.h"

double q_normalisation(double p0, double p1){
    return p1-p0;
}

double p_normalisation(double p0, double p1){
    return p1-p0;
}

double p_weight_function(double p){
    return 1;
}

double q_weight_function(double p){
    return 1;
}

void make_p_bin_quadrature_grids(){

}

void make_q_bin_quadrature_grids(){

}