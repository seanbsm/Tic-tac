/* Spin-independent Malfliet-Tjon potential: short-range repulsive + long-range attractive
 * defined to act only in l=0 states, i.e. it a pure S-wave interaction.
 * This is the MT I-III (spin-singlet/triplet) version with parameters from
 * Tab. 1 PRC 22, 823 (Payne, Friar, Gibson, Afnan)
 *
 *             V_i                 V_i          (p+p')^2 + mu_i^2
 *   ---------------------- =  ----------- ln ----------------------     (in S-waves)
 *   (2pi^2) (q^2 + mu_i^2)      2pi p p'       (p-p')^2 + mu_i^2
 *
 *
 * E_d = 2.23 MeV */

#include "malfliet_tjon.h"

malfliet_tjon::malfliet_tjon(){
}

void malfliet_tjon::update_parameters(double* parameters){
}

void malfliet_tjon::V(double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	Varray[0] = 0;
	Varray[1] = 0;
	Varray[2] = 0;
	Varray[3] = 0;
	Varray[4] = 0;
	Varray[5] = 0;
	
    /* Potential parameters.
     * "R" is for repulsive.
     * "A" is for attractive. */
    double VA_I   =  513.968; // [MeV fm]
    double VA_III =  626.885; // [MeV fm]
    double VR     = 1438.720; // [MeV fm]
    double muA    =    1.550; // [no units]
    double muR    =    3.110; // [no units]

    muA *= hbarc;             // [MeV fm]
    muR *= hbarc;             // [MeV fm]

    /* "num" - numerator
     * "den" - denominator */
    if (J==0){  // MT-I
        double num_A = (qo+qi)*(qo+qi) + muA*muA;
        double num_R = (qo+qi)*(qo+qi) + muR*muR;
        double den_A = (qo-qi)*(qo-qi) + muA*muA;
        double den_R = (qo-qi)*(qo-qi) + muR*muR;

        Varray[0] = ( (VR/hbarc)*std::log(num_R/den_R) - (VA_I/hbarc)*std::log(num_A/den_A) ) / (2*M_PI*qi*qo);
    }
    else if (J==1){  // MT-III
        double num_A = (qo+qi)*(qo+qi) + muA*muA;
        double num_R = (qo+qi)*(qo+qi) + muR*muR;
        double den_A = (qo-qi)*(qo-qi) + muA*muA;
        double den_R = (qo-qi)*(qo-qi) + muR*muR;

        Varray[2] = ( (VR/hbarc)*std::log(num_R/den_R) - (VA_III/hbarc)*std::log(num_A/den_A) ) / (2*M_PI*qi*qo);
    }
}
