#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include "type_defs.h"

/*
 * Here, every single constant we use are listed.
 * This is to have simple overview and prevent mismatch in definition of,
 * for example, proton mass
 */

const floatType deg_to_rad  = M_PI/180.;		     // conversion from degrees to radians	 [rad/deg]
const floatType rad_to_deg  = 180./M_PI;		     // conversion from radians to degress	 [deg/rad]
const floatType divTwo	    = 0.5;				     // factor 0.5
const floatType one		    = 1.;				     // factor 1.0
const floatType two		    = 2.;				     // factor 2.0
const floatType pi		    = M_PI;				     // pi									 [no units]
const floatType Mp		    = 938.272;			     // proton mass  						 [MeV]
const floatType Mn		    = 939.565;			     // neutron mass 						 [MeV]
const floatType MN		    = 0.5*(Mp+Mn);//2*Mp*Mn/(Mp+Mn);	     // Nucleon mass					     [MeV]
const floatType Ed_measured = 2.22457;               // Measured deuteron binding energy     [MeV]
const floatType Md_measured = Mn + Mp - Ed_measured; // Measured deuteron mass               [MeV]
const floatType Lambda	    = 450;				     // cut-off for renormalization of LO 	 [MeV]
const floatType gA 		    = 1.289;			     // axial coupling constant 			 [no units]
const floatType fpi 	    = 92.2;				     // pion decay constant 				 [MeV] 		(used convention 130.41/sqrt(2) = 92.2)
const floatType mpi 	    = 138.039;			     // pion_0 mass 						 [MeV] 		(average of pi+, pi-, and pi0 masses)
const floatType hbarc       = 197.327;               // MeV fm










/* LO_INTERNAL CONSTANTS */
const floatType C1S0	= -0.112927/100.;	// contact term C1S0 for lambda = 450	[MeV]
const floatType C3S1	= -0.087340/100.;	// contact term C3S1 for lambda = 450	[MeV]


/* N2LO_OPT CONSTANTS */
/* LO contacts */
const floatType Ct_1S0pp_n2lo_opt   = -1.51366037000000009e-01; // Ct_1S0pp
const floatType Ct_1S0np_n2lo_opt   = -1.52141088000000008e-01; // Ct_1S0np
const floatType Ct_1S0nn_n2lo_opt   = -1.51764746000000006e-01; // Ct_1S0nn
const floatType Ct_3S1_n2lo_opt     = -1.58434177000000009e-01; // Ct_3S1nn/np/pp (isospin conservation)
  
/* NLO contacts */
const floatType C_1S0_n2lo_opt 	    =  2.40402194400000013e+00; // C_1S0
const floatType	C_3P0_n2lo_opt 	    =  1.26339076300000008e+00; // C_3P0
const floatType	C_1P1_n2lo_opt 	    =  4.17045542000000047e-01; // C_1P1
const floatType	C_3P1_n2lo_opt 	    = -7.82658500000000146e-01; // C_3P1
const floatType	C_3S1_n2lo_opt	    =  9.28384663000000110e-01; // C_3S1
const floatType	C_3S1_3D1_n2lo_opt  =  6.18141418999999970e-01; // C_3S1-3D1
const floatType	C_3P2_n2lo_opt 	    = -6.77808510999999947e-01; // C_3P2

/* Long-range constants */
const floatType c1_n2lo_opt		   	= -9.18639529000000010e-01;
const floatType c3_n2lo_opt		   	= -3.88868749299999994e+00;
const floatType c4_n2lo_opt		   	=  4.31032716100000002e+00;
const floatType gA_n2lo_opt  	    =  1.29000000000000004e+00;



/* IDAHO_N3LO CONSTANTS */
/* LO contacts */
const floatType Ct_1S0pp_idaho_n3lo   	  = -1.45285999999999998e-01; // Ct_1S0pp
const floatType Ct_1S0np_idaho_n3lo   	  = -1.47166999999999992e-01; // Ct_1S0np
const floatType Ct_1S0nn_idaho_n3lo   	  = -1.46284999999999998e-01; // Ct_1S0nn
const floatType Ct_3S1_idaho_n3lo     	  = -1.18972495999999997e-01; // Ct_3S1nn/np/pp (isospin conservation)
  
/* NLO contacts */
const floatType C_1S0_idaho_n3lo 	   	  =  2.37999999999999989e+00; // C_1S0
const floatType	C_3P0_idaho_n3lo 	   	  =  1.48700000000000010e+00; // C_3P0
const floatType	C_1P1_idaho_n3lo 	   	  =  6.56000000000000028e-01; // C_1P1
const floatType	C_3P1_idaho_n3lo 	   	  = -6.30000000000000004e-01; // C_3P1
const floatType	C_3S1_idaho_n3lo	   	  =  7.60000000000000009e-01; // C_3S1
const floatType	C_3S1_3D1_idaho_n3lo   	  =  8.25999999999999956e-01; // C_3S1-3D1
const floatType	C_3P2_idaho_n3lo 	   	  = -5.38000000000000034e-01; // C_3P2
  
/* N3LO contacts */
const floatType	Dh_1S0_idaho_n3lo 	   	  = -2.54499999999999993e+00; // Dh_1S0
const floatType	D_1S0_idaho_n3lo	   	  = -1.60000000000000000e+01; // D_1S0
const floatType	D_3P0_idaho_n3lo	   	  =  2.44999999999999996e-01; // D_3P0
const floatType	D_1P1_idaho_n3lo	   	  =  5.25000000000000000e+00; // D_1P1
const floatType	D_3P1_idaho_n3lo	   	  =  2.35000000000000009e+00; // D_3P1
const floatType	Dh_3S1_idaho_n3lo 	   	  =  7.00000000000000089e+00; // Dh_3S1
const floatType	D_3S1_idaho_n3lo 	   	  =  6.54999999999999982e+00; // D_3S1
const floatType	D_3D1_idaho_n3lo 	   	  = -2.79999999999999982e+00; // D_3D1
const floatType	Dh_3S1_3D1_idaho_n3lo  	  =  2.25000000000000000e+00; // Dh_3S1-3D1
const floatType	D_3S1_3D1_idaho_n3lo   	  =  6.61000000000000032e+00; // D_3S1-3D1
const floatType	D_1D2_idaho_n3lo 	   	  = -1.77000000000000002e+00; // D_1D2
const floatType	D_3D2_idaho_n3lo 	   	  = -1.45999999999999996e+00; // D_3D2
const floatType	D_3P2_idaho_n3lo 	   	  =  2.29499999999999993e+00; // D_3P2
const floatType	D_3P2_3F2_idaho_n3lo   	  = -4.65000000000000024e-01; // D_3P2-3F2
const floatType	D_3D3_idaho_n3lo 	   	  =  5.66000000000000014e+00; // D_3D3

/* Long-range constants */
const floatType c1_idaho_n3lo		   	  = -8.10000000000000053e-01;
const floatType c2_idaho_n3lo		   	  =  2.79999999999999982e+00;
const floatType c3_idaho_n3lo		   	  = -3.20000000000000018e+00;
const floatType c4_idaho_n3lo		   	  =  5.40000000000000036e+00;
const floatType d1_plus_d2_idaho_n3lo 	  =  3.06000000000000005e+00;
const floatType d3_idaho_n3lo         	  = -3.27000000000000002e+00;
const floatType d5_idaho_n3lo         	  =  4.50000000000000011e-01;
const floatType d14_minus_d15_idaho_n3lo  = -5.65000000000000036e+00;
const floatType gA_idaho_n3lo 	  	      =  1.29000000000000004e+00;

#endif // CONSTANTS_H
