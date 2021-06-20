
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <fstream>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "make_pw_symm_states.h"
#include "kinetics.h"
#include "permutation_operators.h"
#include "state_antisymmetrization.h"
#include "faddeev_iterator.h"
#include "General_functions/gauss_legendre.h"
#include "Interactions/potential_model.h"
#include "Triton_states/read_psi.h"

#include "lippmann_schwinger_solver.h"

#include "auxiliary.h"

using namespace std;

void open_file(std::ofstream &file,
			   std::string file_path){
	bool rewrite_file = true;
	if (rewrite_file){
		/* Overwrite file */
		file.open(file_path);
	}
	else{
		/* Append to file */
		file.open(file_path, std::ios_base::app);
	}
}

void lin_interpolate_matrix(double* ref_matrix,
						    double* spln_matrix,
						    double* ref_p_array,
						    double* spln_p_array,
						    int ref_N,
						    int N){

	for (int i = 0; i<N; i++){
		for (int j = 0; j<N; j++){

			double pi = spln_p_array[i];
            double pj = spln_p_array[j];

        	int pi1_index = 0; while ((ref_p_array[pi1_index] < pi) && (pi1_index < ref_N - 1)) pi1_index++; if (pi1_index > 0) pi1_index--;
        	int pj1_index = 0; while ((ref_p_array[pj1_index] < pj) && (pj1_index < ref_N - 1)) pj1_index++; if (pj1_index > 0) pj1_index--;
			
			int pi2_index = pi1_index + 1;
			int pj2_index = pj1_index + 1;
		
			double pi1 = ref_p_array[pi1_index]; double pi2 = ref_p_array[pi2_index];
			double pj1 = ref_p_array[pj1_index]; double pj2 = ref_p_array[pj2_index];

			double M11 = ref_matrix[pi1_index*ref_N + pj1_index];
			double M12 = ref_matrix[pi1_index*ref_N + pj2_index];
			double M21 = ref_matrix[pi2_index*ref_N + pj1_index];
			double M22 = ref_matrix[pi2_index*ref_N + pj2_index];

			double x2x  = pi2 - pi;
			double xx1  = pi  - pi1;
			double x2x1 = pi2 - pi1;
			double y2y  = pj2 - pj;
			double yy1  = pj  - pj1;
			double y2y1 = pj2 - pj1;

			double lin_interpolated_element = (M11*x2x*y2y + M21*xx1*y2y + M12*x2x*yy1 + M22*xx1*yy1)/(x2x1*y2y1);
			
			spln_matrix[i*N + j] = lin_interpolated_element;
		}
	}
}

void triplet_phase_shifts(cfloatType &T11, cfloatType &T12, cfloatType &T22, floatType q, double M, cfloatType* delta_array){
	const floatType radToDeg = 180/pi;
	const cfloatType i{ 0.0, 1.0 };

	cfloatType fac = divTwo*pi*M*q;
	
	/* Blatt-Biedenharn (BB) convention */
	cfloatType twoEpsilonJ_BB = atan(two*T12/(T11-T22));	// mixing parameter
	cfloatType delta_plus_BB  = -divTwo*i*log(one - i*fac*(T11+T22) + i*fac*(two*T12)/sin(twoEpsilonJ_BB));
	cfloatType delta_minus_BB = -divTwo*i*log(one - i*fac*(T11+T22) - i*fac*(two*T12)/sin(twoEpsilonJ_BB));
	
	/* Stapp convention (bar-phase shifts) in terms of Blatt-Biedenharn convention (New formula, no kink) */
	floatType cos2eps = cos(divTwo*twoEpsilonJ_BB.real())*cos(divTwo*twoEpsilonJ_BB.real());
	floatType cos_2delta_plus  = cos(2.*delta_plus_BB.real());
	floatType sin_2delta_plus  = sin(2.*delta_plus_BB.real());
	floatType cos_2delta_minus = cos(2.*delta_minus_BB.real());
	floatType sin_2delta_minus = sin(2.*delta_minus_BB.real());
	
	floatType aR, aI, tmp;
	
	aR = cos2eps*cos_2delta_minus + (1.-cos2eps)*cos_2delta_plus;
	aI = cos2eps*sin_2delta_minus + (1.-cos2eps)*sin_2delta_plus;
	
	cfloatType delta_minus = 0.5*atan2(aI, aR) * radToDeg;
	
	aR = cos2eps*cos_2delta_plus + (1.-cos2eps)*cos_2delta_minus;
	aI = cos2eps*sin_2delta_plus + (1.-cos2eps)*sin_2delta_minus;
	
	cfloatType delta_plus = 0.5*atan2(aI, aR) * radToDeg;
	
	tmp = 0.5*sin(twoEpsilonJ_BB.real());
	aR  = tmp*(cos_2delta_minus - cos_2delta_plus);
	aI  = tmp*(sin_2delta_minus - sin_2delta_plus);
	tmp = (delta_plus.real() + delta_minus.real())/radToDeg;
	
	cfloatType epsilon = 0.5*asin(aI*cos(tmp) - aR*sin(tmp)) * radToDeg;
	
	delta_array[0] = delta_minus;
	delta_array[1] = delta_plus;
	delta_array[2] = epsilon;
}

void update_potential_matrix(double* V_array,
							 double* p_array,
							 double E,
							 int Nk,
							 int L, int Lp, int S, int J, int T,
							 bool coupled,
							 potential_model* pot_ptr){
	int Nk1 = Nk + 1;
	int Nk2 = 2*Nk1;

	bool E_positive = (E>=0);
	double p = 0;
	if (E_positive){
		p = sqrt(E*MN);
	}
	
	/* Construct 2N potential matrix <k|v|k_p> */
	double V_elements_col [6];
	double V_elements_row [6];
	double k=0, k_in=0, k_p=0, k_out=0;
	for (int i=0; i<Nk; i++){

        k  = p_array[i];//*hbarc;

        if(E_positive){
			pot_ptr->V(k, p, coupled, S, J, T, V_elements_col);
			pot_ptr->V(p, k, coupled, S, J, T, V_elements_row);
		}
		else{
			for(int i=0;i<6;i++){
				V_elements_col[i]=0;
				V_elements_row[i]=0;
			}
		}

		if (coupled){
			/* Set end-column <ki|v|p> */
			V_array[ i	   *Nk2 +  Nk] 		= extract_potential_element_from_array(J-1, J-1, J, S, coupled, V_elements_col);
			V_array[ i	   *Nk2 + (Nk+Nk1)] = extract_potential_element_from_array(J-1, J+1, J, S, coupled, V_elements_col);
			V_array[(i+Nk1)*Nk2 +  Nk] 		= extract_potential_element_from_array(J+1, J-1, J, S, coupled, V_elements_col);
			V_array[(i+Nk1)*Nk2 + (Nk+Nk1)] = extract_potential_element_from_array(J+1, J+1, J, S, coupled, V_elements_col);
			/* Set end-row <p'|v|ki> */
			V_array[ Nk		*Nk2 +  i] 		= extract_potential_element_from_array(J-1, J-1, J, S, coupled, V_elements_row);
			V_array[ Nk		*Nk2 + (i+Nk1)] = extract_potential_element_from_array(J-1, J+1, J, S, coupled, V_elements_row);
			V_array[(Nk+Nk1)*Nk2 +  i] 		= extract_potential_element_from_array(J+1, J-1, J, S, coupled, V_elements_row);
			V_array[(Nk+Nk1)*Nk2 + (i+Nk1)] = extract_potential_element_from_array(J+1, J+1, J, S, coupled, V_elements_row);
		}
		else{
			/* Set end-column <ki|v|p> */
			V_array[i*Nk1 + Nk] = extract_potential_element_from_array(L, Lp, J, S, coupled, V_elements_col);
			/* Set end-row <p'|v|ki> */
			V_array[Nk*Nk1 + i] = extract_potential_element_from_array(L, Lp, J, S, coupled, V_elements_row);
		}
	}

	if(E_positive){pot_ptr->V(p, p, coupled, S, J, T, V_elements_row);}
	else{ for(int i=0;i<6;i++){ V_elements_row[i]=0; } }
	
	if (coupled){
		V_array[ Nk	    *Nk2 +  Nk] 	 = extract_potential_element_from_array(J-1, J-1, J, S, coupled, V_elements_row);
		V_array[ Nk	    *Nk2 + (Nk+Nk1)] = extract_potential_element_from_array(J-1, J+1, J, S, coupled, V_elements_row);
		V_array[(Nk+Nk1)*Nk2 +  Nk] 	 = extract_potential_element_from_array(J+1, J-1, J, S, coupled, V_elements_row);
		V_array[(Nk+Nk1)*Nk2 + (Nk+Nk1)] = extract_potential_element_from_array(J+1, J+1, J, S, coupled, V_elements_row);
	}
	else{
		/* Set last element <p'|v|p> */
		V_array[Nk1*Nk1 - 1] = extract_potential_element_from_array(L, Lp, J, S, coupled, V_elements_row);
	}
}

void store_array(double *arrays,
				 int num_rows,
				 int num_cols,
				 std::string folder_name){
	
	std::string file_path  = folder_name + ".csv";
	
	/* Open file*/
	std::ofstream result_file;
	open_file(result_file, file_path);
	
	/* Fixes formatting of stored numbers */
	result_file << std::fixed
				<< std::showpos
				<< std::setprecision(8);
	
	/* Append cross-sections */
	for (int i=0; i<num_rows; i++){
		for (int j=0; j<num_cols; j++){
			/* Append vector element */
			result_file << arrays[i*num_cols + j] << ",";
		}
		result_file << "\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}


int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* Start main body of code here */

	/* Tritium bound-state quantum numbers */
	int two_J_3N  	 = 1;
    int two_T_3N  	 = 1;
    int parity_3N 	 = 1;

	/* PWE truncation */
	/* Maximum (max) and minimum (min) values for J_2N and J_1N */
    int J_2N_min 	 = 0;	// The LS-solver will fail if this is not zero - I haven't taken this into account in my indexing
    int J_2N_max 	 = 1;

	/* Quadrature 3N momenta */
	int Np	   		 = 32;
	int Nq	   		 = 30;
	int Nx 			 = 20;
	int Nalpha 		 =  0;
	double* p_array  = NULL;
	double* q_array  = NULL;
	double* wp_array = NULL;
	double* wq_array = NULL;
    double* x_array  = new double [Nx];
    double* wx_array = new double [Nx];
	calculate_angular_quadrature_grids(x_array, wx_array, Nx);

	/* Momentum-representation of 3N state at quadrature nodes in p_array and q_array */
	double* state_3N_symm_array = NULL;
	double* state_3N_asym_array = NULL;

	/* Permutation operator */
	double* P123_array  = NULL;

	/* Anti-symmetrization operator */
	double* A123  = NULL;

	/* Quantum numbers of partial-wave expansion in state_3N_array */
    int* L_2N     = NULL;     // pair angular momentum
    int* S_2N     = NULL;     // pair total spin
    int* J_2N     = NULL;     // pair total angular momentum
    int* T_2N     = NULL;     // pair total isospin
    int* l_3N     = NULL;     // three-nucleon angular momentum (?)
    int* two_j_3N = NULL; 	  // three-nucleon total angular momentum x2 (?)

	/* Quantum numbers of partial-wave expansion in P123 */
    int* L_2N_P123     = NULL;     // pair angular momentum
    int* S_2N_P123     = NULL;     // pair total spin
    int* J_2N_P123     = NULL;     // pair total angular momentum
    int* T_2N_P123     = NULL;     // pair total isospin
    int* l_3N_P123     = NULL;     // three-nucleon angular momentum (?)
    int* two_j_3N_P123 = NULL; 	   // three-nucleon total angular momentum x2 (?)

	/* Tells the program to read pre-calculated antisymmetric triton states.
	 * Handy for small tests since the P123-file can be huge */
	bool use_premade_symmetric_states 	  = false;
	bool use_premade_antisymmetric_states = false;

	//potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("LO_internal", "np");
	//potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("LO_internal", "nn");
	//potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("N2LOopt", "np");
	//potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("N2LOopt", "nn");
	potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_N3LO", "np");
	potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_N3LO", "nn");
	//potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("Idaho_EM500", "np");
	//potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("Idaho_EM500", "nn");

	if (use_premade_symmetric_states){
		cout << "Reading states from file" << endl;
		get_all_states(&state_3N_symm_array, &state_3N_asym_array, Np, &p_array, &wp_array, Nq, &q_array, &wq_array, Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);
	}
	else{
		cout << "Constructing 3N partial-wave basis" << endl;
		construct_symmetric_pw_states(two_J_3N, two_T_3N, parity_3N,
									  J_2N_min, J_2N_max,
									  Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);
		p_array  = new double [Np];
		wp_array = new double [Np];
		q_array  = new double [Nq];
		wq_array = new double [Nq];
		
		bool use_trigonometric_distribution = false;
        if (use_trigonometric_distribution){
            cout << "Constructing p mesh" << endl;
            gauss(p_array, wp_array, Np); rangeChange_0_inf(p_array, wp_array, 1000., Np);
	        cout << "Constructing q mesh" << endl;
	        gauss(q_array, wq_array, Nq); rangeChange_0_inf(q_array, wq_array, 1000., Nq);
        }
        else{
            double min = 0; double max=5;
            cout << "Constructing p mesh" << endl;
            calc_gauss_points (p_array, wp_array, min, max, Np);
	        cout << "Constructing q mesh" << endl;
	        calc_gauss_points (q_array, wq_array, min, max, Nq);

			/* Convert from fm^-1 to MeV */
			for (int i = 0; i<Np; i++){
				p_array[i]  *= hbarc;
				wp_array[i] *= hbarc;
			}
			for (int i = 0; i<Nq; i++){
				q_array[i]  *= hbarc;
				wq_array[i] *= hbarc;
			}
        }
	}

	if (use_premade_antisymmetric_states == true){
		cout << "Reading P123 from file and calculating A123" << endl;
		calculate_antisymmetrization_operator(Np, Nq, Nalpha, &A123, q_array, p_array);

		cout << "Calculating psi_asym from psi_symm using A123" << endl;
		antisymmetrize_state(state_3N_symm_array, state_3N_asym_array, A123, Np, Nq, Nalpha);
	}
	else{
		cout << "Read P123 dimensions from h5-file" << endl;
		int Nalpha_P123 = 0;
		int Np_P123 = 0;
		int Nq_P123 = 0;
		get_h5_P123_dimensions(Nalpha_P123, Np_P123, Nq_P123);
    	P123_array = new double [Np_P123 * Nq_P123 * Nalpha_P123 * Np_P123 * Nq_P123 * Nalpha_P123];
		
		int* L_2N_P123     = new int [Nalpha_P123];     // pair angular momentum
    	int* S_2N_P123     = new int [Nalpha_P123];     // pair total spin
    	int* J_2N_P123     = new int [Nalpha_P123];     // pair total angular momentum
    	int* T_2N_P123     = new int [Nalpha_P123];     // pair total isospin
    	int* l_3N_P123     = new int [Nalpha_P123];     // three-nucleon angular momentum (?)
    	int* two_j_3N_P123 = new int [Nalpha_P123]; 	// three-nucleon total angular momentum x2 (?)

		if (Np!=Np_P123){
			raise_error("P123 Np-dimension does not match grid setup");
		}
		if (Nq!=Nq_P123){
			raise_error("P123 Nq-dimension does not match grid setup");
		}
		if (Nalpha>Nalpha_P123){
			raise_error("P123 Nalpha-dimension does not match grid setup");
		}

		cout << "Reading P123 from file" << endl;
		read_P123_h5_data_file(P123_array,
							   Nq, q_array,
							   Np, p_array,
							   Nalpha_P123, L_2N_P123, S_2N_P123, J_2N_P123, T_2N_P123, l_3N_P123, two_j_3N_P123);

		///* Convert from fm^-1 to MeV */
		//for (int i = 0; i<Np; i++){
		//	p_array[i]  *= hbarc;
		//	wp_array[i] *= hbarc;
		//}
		//for (int i = 0; i<Nq; i++){
		//	q_array[i]  *= hbarc;
		//	wq_array[i] *= hbarc;
		//}

		/* Check that the P123 PW states have the same ordering as the local PW states */
		//for (int i=0; i<Nalpha; i++){
		//	cout << (L_2N_P123[i]      == L_2N[i]    ) << endl;
    	//	cout << (S_2N_P123[i]      == S_2N[i]    ) << endl;    
    	//	cout << (J_2N_P123[i]      == J_2N[i]    ) << endl;   
    	//	cout << (T_2N_P123[i]      == T_2N[i]    ) << endl;    
    	//	cout << (l_3N_P123[i]      == l_3N[i]    ) << endl;   
    	//	cout << (two_j_3N_P123[i]  == two_j_3N[i]) << endl;
		//}


		/* Convert from fm^-1 to MeV */
		//for (int i = 0; i<Np; i++){
		//	p_array[i]  *= hbarc;
		//	wp_array[i] *= hbarc;
		//}
		//for (int i = 0; i<Nq; i++){
		//	q_array[i]  *= hbarc;
		//	wq_array[i] *= hbarc;
		//}

		/*for (int idx_p = 0; idx_p<Np; idx_p++){
			p_array[idx_p] /= hbarc;
		}
		for (int idx_q = 0; idx_q<Nq; idx_q++){
			q_array[idx_q] /= hbarc;
		}*/

		/* Test diagonalization routine */
		/*cout << "Start test:" << endl;
		double* A = new double [2*2];
		A[0] = 2;
		A[1] = 1;
		A[2] = 1;
		A[3] = 2;
		double* v = new double [2*2];
		double* w = new double [2];
		find_eigenvalues(A, w, v, 2);
		for (int i = 0; i<2; i++){
			cout << w[i] << " ";
			for (int j = 0; j<2; j++){
				cout << v[j*2 + i] << " ";
			}
			cout << endl;
		}
		cout << "End test:" << endl;*/
		
		/* Test momentum grid */
		/*int N=100;
		double* p_par = new double [N]; double* p_mom = new double [N];
		double* w_par = new double [N];
		double min = 0; double max=6.5;
		calc_gauss_points (p_par, w_par, min, max, N);
		for (int i = 0; i<N; i++){
			//printf("%.5e %.5e \n", w_par[i], p_par[i]);
			p_mom[i] = p_par[i]*hbarc;
		}*/

		/* Test NN matrix elements */
		/*
		double Vnn [6]; double Vnp [6]; double V [6];
		double p0 = 1e-15;
		int nul = 0;
		int one = 1;
		pot_ptr_np->V(p0,p0,false,nul,nul,one,V);
		pot_ptr_nn->V(p0, p0, false, nul, nul, one, Vnn);
        pot_ptr_np->V(p0, p0, false, nul, nul, one, Vnp);
	    for (int i=0; i<6; i++){
	    	V[i] = (1./3)*Vnp[i] + (2./3)*Vnn[i];
	    }
		std::cout << Vnp[0] << " " << hbarc*hbarc*Vnp[0] << std::endl;
		return 0;
		*/

		/* Advanced test NN matrix elements */
		/*double V [6];
		int S = 0;
		int L = 0; int Lp = 0;
		int J = 0;
		int T = 1;
		double largest;
	    for (int i = 0; i<N; i++){
			for (int j = 0; j<N; j++){
				pot_ptr_np->V(p_mom[i], p_mom[j],false,S,J,T,V);
				//if (abs(V1S0[i*N + j])>1e-7){
				//	printf("%.5e %.5e %.5e %.5e %.5e\n", p_par[i], p_par[j], hbarc*MN*V[0]*M_PI/2, V1S0[i*N + j], abs(hbarc*MN*V[0]*M_PI/2-V1S0[i*N + j]));
				//}
				//if (abs(p_mom[i]-p_mom[j])>5e2){
				//	printf("%.5e %.5e %.5e %.5e %.5e\n", p_par[i], p_par[j], hbarc*MN*V[0]*M_PI/2, V1S0[i*N + j], abs(p_mom[i]-p_mom[j]));
				//}
				//if (V[0]==0 && abs(V1S0[i*N + j])>largest){
				//	largest = V1S0[i*N + j];
				//}
				
				//if (V[0]==0 && p_par[i]<5 && p_par[j]<5){
				//	cout << "hey!" << endl;
				//}
				
				printf("%.5e %.5e %.5e\n", p_par[i], p_par[j], hbarc*MN*V[0]*M_PI/2);
				V[0] = 0;
			}
		}
		cout << largest << endl;
		return 0;*/

		/* Test T-matrix elements */
		/*int N=30;
		double* p_par = new double [N]; double* p_mom = new double [N];
		double* w_par = new double [N]; double* w_mom = new double [N];
		double min = 0; double max=5;
		//gauss(p_mom, w_mom, N);
		//rangeChange_0_inf(p_mom, w_mom, 1000., N);
		//rangeChange_0_inf(p_mom, w_mom, max*hbarc, N);
		//updateRange_a_b(p_mom, w_mom, min*hbarc, max*hbarc, N);
		calc_gauss_points (p_par, w_par, min, max, N);
		for (int i = 0; i<N; i++){
			p_mom[i] = p_par[i]*hbarc;
			w_mom[i] = w_par[i]*hbarc;
			//p_par[i] = p_mom[i]/hbarc;
			//w_par[i] = w_mom[i]/hbarc;
		}

		double* V_unco_array = new double [N*N   * 2*J_2N_max];
		double* V_coup_array = new double [N*N*4 *   J_2N_max];
		calculate_potential_matrices_array(V_unco_array,
                                           V_coup_array,
                                           N, p_mom, w_par,
                                           Nalpha, L_2N, S_2N, J_2N, T_2N,
                                           pot_ptr_nn,
                                           pot_ptr_np);

		for (int i = 0; i<N; i++){
			for (int j = 0; j<N; j++){
				V_unco_array[i*N+j] *= hbarc*MN*M_PI/2;
			}
		}

		//int N_p=100;
		//double* p = new double [N_p];
		//double* w = new double [N_p];
		//min = 0; max=6.5;
		//calc_gauss_points (p, w, min, max, N_p);
		double* V_1S0_interpolated = new double [N*N];
		
		double extern_to_local_conversion = 2/(hbarc*MN*M_PI);
		double local_to_extern_conversion = hbarc*MN*M_PI/2;
		for (int i = 0; i<N; i++){
			//printf("%.8f %.8f \n", p_mom[i], w_mom[i]);
			for (int j = 0; j<N; j++){
				int S=0;
				int J=0;
				int T=1; double V [6];
				pot_ptr_np->V(p_mom[i], p_mom[j],false,S,J,T,V);

				V_1S0_interpolated[i*N + j] = V[0];
				V[0] = 0;
			}
		}

		double T_lab = 10;
		double q_com = sqrt( Mn*Mn*T_lab*(T_lab + 2*Mp) / ( (Mp+Mn)*(Mp+Mn) + 2*T_lab*Mn ) );
		double Z1 = -7.5; double Z2 = -8;//q_com*q_com/MN;
		double q_fm = 0;//0.00854292;
		double q_MeV = q_fm*hbarc;
		//double E_LS1 = Z1*MN/(hbarc*hbarc)  - 0.75*q_fm*q_fm; double E_LS2 = Z2*MN/(hbarc*hbarc)  - 0.75*q_fm*q_fm;
		double E_LS1 = Z1  - 0.75*q_MeV*q_MeV/MN; double E_LS2 = Z2  - 0.75*q_MeV*q_MeV/MN;

		double* V_array = new double [(N+1)*(N+1)];
		for (int i = 0; i<N; i++){
			for (int j = 0; j<N; j++){
				V_array[i*(N+1)+j] = V_1S0_interpolated[i*N+j];
			}
		}
		update_potential_matrix(V_array,
								p_mom,
								E_LS2,
								N,
								0,0,0,0,1,
								false,
								pot_ptr_np);
		
		cfloatType* t_unco_array1_complex = new cfloatType [(N+1)*(N+1)];
		cfloatType* t_unco_array2_complex = new cfloatType [(N+1)*(N+1)];
		calculate_t_element(V_array,
                            t_unco_array1_complex,
                            false,
		                	E_LS1, MN,
		                	N, p_mom, w_mom,
							0, 0);
		calculate_t_element(V_array,
                            t_unco_array2_complex,
                            false,
		                	E_LS2, MN,
		                	N, p_mom, w_mom,
							0, 0);*/
		
		//double* t_unco_array1 = new double [(N+1)*(N+1)];
		//double* t_unco_array2 = new double [(N+1)*(N+1)];
		//for (int i=0; i<(N+1)*(N+1); i++){
		//	t_unco_array1[i] = t_unco_array1_complex[i].real();
		//	t_unco_array2[i] = t_unco_array2_complex[i].real();
		//}

		//floatType b = 2*(Mp - q_com*q_com/Mn);
		//floatType c = -q_com*q_com*(Mp+Mn)*(Mp+Mn)/(Mn*Mn);
		//T_lab = 0.5*(-b + sqrt(b*b - 4*c));
		//printf("%.5f \n", T_lab);

		//if (E_LS2>=0){
		//	const floatType divTwo	= 0.5;				// factor 0.5
		//	const floatType one		= 1.;				// factor 1.0
		//	const floatType two		= 2.;				// factor 2.0
		//	const floatType radToDeg = 180/pi;
		//	const cfloatType I {0.0, 1.0};
		//	double q_on_shell = sqrt(E_LS2*MN);
		//	cfloatType T_on_shell = t_unco_array2_complex[(N+1)*(N+1)-1];
		//	cfloatType Z = -MN*q_on_shell*I*T_on_shell*M_PI + one;
		//	cfloatType delta = (-divTwo*I*log(Z)*radToDeg).real();
		//	//printf(delta, "\n");
		//	printf("%.5f %.5f\n", delta.real(), delta.imag());
		//}

		//int spln_N=30;
		//double* spln_p_array_fm = new double [spln_N]; double* spln_p_array_MeV = new double [spln_N];
		//double* spln_w_array_fm = new double [spln_N]; double* spln_w_array_MeV = new double [spln_N];
		//min = 0; max=5;
		//calc_gauss_points (spln_p_array_fm, spln_w_array_fm, min, max, spln_N);
		//for (int i = 0; i<N; i++){
		//	spln_p_array_MeV[i] = spln_p_array_fm[i]*hbarc;
		//	spln_w_array_MeV[i] = spln_w_array_fm[i]*hbarc;
		//}
		//double* t_mat = new double [spln_N*spln_N]; 
		//lin_interpolate_matrix(t_unco_array2,
		//				       t_mat,
		//				       p_mom,
		//				       spln_p_array_MeV,
		//				       N,
		//				       spln_N);
		
		// PRINTS SPLINED T-MATRIX
		//for (int i = 0; i<1; i++){
		//	for (int j = 0; j<spln_N; j++){
		//		printf("%.5e %.5e %.5e\n", spln_p_array_fm[i], spln_p_array_fm[j], t_mat[i*spln_N+j]*local_to_extern_conversion);
		//	}
		//}

		// PRINTS EXACT T-MATRIX
		//for (int i = 0; i<1; i++){
		//	for (int j = 0; j<N; j++){
		//		printf("%.5e %.5e %.5e\n", p_par[i], p_par[j], t_unco_array2[i*(N+1)+j]*local_to_extern_conversion);
		//	}
		//}

		// PRINTS TWO T-MATRICES (EXACT)
		//for (int i = 0; i<1; i++){
		//	for (int j = 0; j<N; j++){
		//		printf("%.5e %.5e %.5e %.5e\n", p_par[i], p_par[j], t_unco_array1[i*(N+1)+j]*local_to_extern_conversion, t_unco_array2[i*(N+1)+j]*local_to_extern_conversion);
		//	}
		//}

		//double* arrays = new double [N*N*3];
		//for (int i=0; i<N; i++){
		//	for (int j=0; j<N; j++){
		//		arrays[(i*N+j)*3 + 0] = p_par[i];
		//		arrays[(i*N+j)*3 + 1] = p_par[j];
		//		//arrays[(i*N+j)*3 + 2] = V_array[i*(N+1)+j]*local_to_extern_conversion;//
		//		arrays[(i*N+j)*3 + 2] = t_unco_array1[i*(N+1)+j]*hbarc*MN;
		//	}
		//}
		//store_array(arrays, N*N, 3, "t_elements");
		////store_array(arrays, N*N, 3, "V_elements");
		//return 0;

		/* CHECK POTENTIAL ARRAYS */
		/*int N=100;
		double* p_par = new double [N]; double* p_mom = new double [N];
		double* w_par = new double [N]; double* w_mom = new double [N];
		double min = 0; double max=6.5;
		calc_gauss_points (p_par, w_par, min, max, N);
		for (int i = 0; i<N; i++){
			p_mom[i] = p_par[i]*hbarc;
			w_mom[i] = w_par[i]*hbarc;
			//p_par[i] = p_mom[i]/hbarc;
			//w_par[i] = w_mom[i]/hbarc;
		}
		int Np1 = N+1;
    	double* V_unco_array = new double [Np1*Np1   * 2*J_2N_max];
    	double* V_coup_array = new double [Np1*Np1*4 *   J_2N_max];

    	calculate_potential_matrices_array(V_unco_array,
    	                                   V_coup_array,
    	                                   N, p_mom, w_mom,
    	                                   Nalpha, L_2N, S_2N, J_2N, T_2N,
    	                                   pot_ptr_nn,
    	                                   pot_ptr_np);
		double* arrays = new double [N*N*3];
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				arrays[(i*N+j)*3 + 0] = p_par[i];
				arrays[(i*N+j)*3 + 1] = p_par[j];
				arrays[(i*N+j)*3 + 2] = V_unco_array[i*(N+1)+j]*hbarc*MN*M_PI/2;
			}
		}
		store_array(arrays, N*N, 3, "V_elements");
		return 0;*/

		/* CHECK COUPLED T-MATRIX */
		/*int N=4;
		double* p_mom = new double [N];
		double* w_mom = new double [N];
		gauss(p_mom, w_mom, N);
		rangeChange_0_inf(p_mom, w_mom, 1000., N);
		int Np1 = N+1;
    	int V_unco_array_size = Np1*Np1   * 2*(J_2N_max+1);
    	int V_coup_array_size = Np1*Np1*4 *   (J_2N_max+1);
    	double* V_unco_array = new double [V_unco_array_size];
    	double* V_coup_array = new double [V_coup_array_size];
		printf("Potential arrays start \n");
		calculate_potential_matrices_array(V_unco_array,
                                           V_coup_array,
                                           N, p_mom, w_mom,
                                           Nalpha, L_2N, S_2N, J_2N, T_2N,
                                           pot_ptr_nn,
                                           pot_ptr_np);
		printf("Potential arrays end \n");
		double T_lab = 10;
		double q_com = sqrt( Mn*Mn*T_lab*(T_lab + 2*Mp) / ( (Mp+Mn)*(Mp+Mn) + 2*T_lab*Mn ) );
		double E = q_com*q_com/MN;
		printf("Potential arrays update start \n");
		update_potential_matrix(V_coup_array,
								p_mom,
								E,
								N,
								0,2,1,1,0,	//int L, int Lp, int S, int J, int T
								true,
								pot_ptr_np);
		printf("Potential arrays update end \n");

		cfloatType* t_coup_array_complex = new cfloatType [4*(N+1)*(N+1)];
		printf("T-matrix arrays start \n");
		calculate_t_element(V_coup_array,
                            t_coup_array_complex,
                            true,
		                	E, MN,
		                	N, p_mom, w_mom);
		printf("T-matrix arrays end \n");

		printf("%.5e \n", q_com);
		for (int i = 0; i<N+1; i++){
			for (int j = 0; j<N+1; j++){
				printf("%.5e %.5e %.5e %.5e\n", p_mom[i], p_mom[j], t_coup_array_complex[i*2*(N+1)+j].real(), t_coup_array_complex[i*2*(N+1)+j].imag());
				printf("%.5e %.5e %.5e %.5e\n", p_mom[i], p_mom[j], t_coup_array_complex[i*2*(N+1)+j+N+1].real(), t_coup_array_complex[i*2*(N+1)+j+N+1].imag());
				printf("%.5e %.5e %.5e %.5e\n", p_mom[i], p_mom[j], t_coup_array_complex[(i+N+1)*2*(N+1)+j].real(), t_coup_array_complex[(i+N+1)*2*(N+1)+j].imag());
				printf("%.5e %.5e %.5e %.5e\n", p_mom[i], p_mom[j], t_coup_array_complex[(i+N+1)*2*(N+1)+j+N+1].real(), t_coup_array_complex[(i+N+1)*2*(N+1)+j+N+1].imag());
				printf("\n");
			}
		}

		if (true){
			cfloatType delta_array [3];
			triplet_phase_shifts(t_coup_array_complex[N*2*(N+1)+N], t_coup_array_complex[N*2*(N+1)+N+N+1], t_coup_array_complex[(N+N+1)*2*(N+1)+N+N+1], q_com, MN, delta_array );
			printf("%.5e %.5e %.5e\n", delta_array[0].real(), delta_array[1].real(), delta_array[2].real());
		}
		return 0;*/



		cout << "Starting Faddeev Iterator" << endl;
		calculate_faddeev_convergence(state_3N_asym_array,
									  P123_array, Nalpha_P123, Np_P123, Nq_P123,
                     				  Np, p_array, wp_array,
                     				  Nq, q_array, wq_array,
                     				  Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N,
                     				  two_T_3N, two_J_3N, parity_3N, J_2N_max,
                     				  pot_ptr_nn,
                     				  pot_ptr_np);
	}
	
	/* This function will calculate the 3N c.m. kinetic energy T */
	cout << "Calculating kinetic energy" << endl;
	double kinetic_energy = calculate_3N_kinetic_energy(state_3N_asym_array, state_3N_asym_array, Np, p_array, wp_array, Nq, q_array, wq_array, Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N);
	
	cout << "Calculating potential energy" << endl;
	double potential_energy = calculate_3N_potential_energy(state_3N_asym_array, state_3N_asym_array, Np, p_array, wp_array, Nq, q_array, wq_array, Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N, pot_ptr_np, pot_ptr_nn);

	cout << "T: " << kinetic_energy << endl;
	cout << "V: " << potential_energy << endl;
	cout << "H: " << kinetic_energy + potential_energy << endl;

	/* End main body of code here */

	auto program_end= chrono::system_clock::now();
	

	chrono::duration<double> total_time = program_end - program_start;
	cout << "Total run-time:                         " << total_time.count() << endl;

	std::cout << "END OF RUN" << std::endl;
	
    return 0;
}

