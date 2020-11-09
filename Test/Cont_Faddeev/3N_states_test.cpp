
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

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

int main(int argc, char* argv[]){

	auto program_start = chrono::system_clock::now();
	
	/* Start main body of code here */

	/* Tritium bound-state quantum numbers */
	int two_J_3N  	 = 1;
    int two_T_3N  	 = 1;
    int parity_3N 	 = 1;

	/* PWE truncation */
	/* Maximum (max) and minimum (min) values for J_2N and J_1N */
    int two_J_1N_min = 1;
    int two_J_1N_max = 4;
    int J_2N_min 	 = 0;	// The LS-solver will fail if this is not zero - I haven't taken this into account in my indexing
    int J_2N_max 	 = 3;

	/* Quadrature 3N momenta */
	int Np	   		 = 32;
	int Nq	   		 = 10;
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

	/* Tells the program to read pre-calculated antisymmetric triton states.
	 * Handy for small tests since the P123-file can be huge */
	bool use_premade_symmetric_states 	  = false;
	bool use_premade_antisymmetric_states = false;

	//potential_model* pot_ptr_np = potential_model::fetch_potential_ptr("LO_internal", "np");
	//potential_model* pot_ptr_nn = potential_model::fetch_potential_ptr("LO_internal", "nn");
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
									  two_J_1N_min, two_J_1N_max, J_2N_min, J_2N_max,
									  Nalpha, &L_2N, &S_2N, &J_2N, &T_2N, &l_3N, &two_j_3N);
		cout << "Constructing p mesh" << endl;
		p_array  = new double [Np];
		wp_array = new double [Np];
		gauss(p_array, wp_array, Np); rangeChange_0_inf(p_array, wp_array, 1000., Np);
		cout << "Constructing q mesh" << endl;
		q_array  = new double [Nq];
		wq_array = new double [Nq];
		gauss(q_array, wq_array, Nq); rangeChange_0_inf(q_array, wq_array, 1000., Nq);
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

		cout << "Reading P123 from file" << endl;
		read_P123_h5_data_file(P123_array, Nq, q_array, Np, p_array);

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
		double* w_par = new double [N];
		double min = 0; double max=5;
		calc_gauss_points (p_par, w_par, min, max, N);
		for (int i = 0; i<N; i++){
			//cout << p_par[i] << endl;
			p_mom[i] = p_par[i]*hbarc;
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

		int N_p=100;
		double* p = new double [N_p];
		double* w = new double [N_p];
		min = 0; max=6.5;
		calc_gauss_points (p, w, min, max, N_p);
		double* V_1S0_interpolated = new double [N*N];
		
		double extern_to_local_conversion = 2/(hbarc*MN*M_PI);
		double local_to_extern_conversion = hbarc*MN*M_PI/2;
		for (int i = 0; i<N; i++){
			for (int j = 0; j<N; j++){

				double pi = p_par[i];
                double pj = p_par[j];

            	int pi1_index = 0; while ((p[pi1_index] < pi) && (pi1_index < N_p - 1)) pi1_index++; if (pi1_index > 0) pi1_index--;
            	int pj1_index = 0; while ((p[pj1_index] < pj) && (pj1_index < N_p - 1)) pj1_index++; if (pj1_index > 0) pj1_index--;
				
				int pi2_index = pi1_index + 1;
				int pj2_index = pj1_index + 1;
			
				double pi1 = p[pi1_index]; double pi2 = p[pi2_index];
				double pj1 = p[pj1_index]; double pj2 = p[pj2_index];

				double V11 = V1S0[pi1_index*N_p + pj1_index];
				double V12 = V1S0[pi1_index*N_p + pj2_index];
				double V21 = V1S0[pi2_index*N_p + pj1_index];
				double V22 = V1S0[pi2_index*N_p + pj2_index];

				double x2x  = pi2 - pi;
				double xx1  = pi  - pi1;
				double x2x1 = pi2 - pi1;
				double y2y  = pj2 - pj;
				double yy1  = pj  - pj1;
				double y2y1 = pj2 - pj1;

				double lin_interpolate = (V11*x2x*y2y + V21*xx1*y2y + V12*x2x*yy1 + V22*xx1*yy1)/(x2x1*y2y1);
				
				
				int S=0;
				int J=0;
				int T=0; double V [6];
				pot_ptr_np->V(p_mom[i], p_mom[j],false,S,J,T,V);

				//V_1S0_interpolated[i*N + j] = V[0];
				V_1S0_interpolated[i*N + j] = lin_interpolate * extern_to_local_conversion;

				//if (pi>5 or pj>5){
				//	V_1S0_interpolated[i*N + j] = lin_interpolate * extern_to_local_conversion;
				//}
				//else{
				//	V_1S0_interpolated[i*N + j] = V[0] ;//* local_to_extern_conversion;
				//}

				//V_1S0_interpolated[i*N + j] = lin_interpolate * 2/(hbarc*MN*M_PI);

				//printf("%.5e %.5e        %.5e %.5e %.5e %.5e       %.5e %.5e %.5e\n", p_par[i], p_par[j], V11, V12, V21, V22, lin_interpolate, V_1S0_interpolated[i*N + j], V[0] * local_to_extern_conversion);
				//printf("%.5e %.5e %.5e\n", p_par[i], p_par[j], V_1S0_interpolated[i*N + j]);
				V[0] = 0;
			}
		}

		double* t_unco_array = new double [  (N+1)*(N+1)];
		double Z = -8;
		double q = 0.00854292*hbarc;
		double E_LS = Z - 0.75*q*q/MN;
		//double E_LS = Z*MN/(hbarc*hbarc) - 0.75*q*q;
		//double t_element = calculate_t_element(&(V_unco_array[0]),
        //                                       t_unco_array,
        //                                       0, 0, 0, 0, 1,
		//                			           E_LS, MN,
		//                			           N, p_mom, w_par,
		//                			           0, 0,
		//                			           pot_ptr_nn, pot_ptr_np);
		double t_element = calculate_t_element(V_1S0_interpolated,
                                               t_unco_array,
                                               0, 0, 0, 0, 1,
		                			           E_LS, MN,
		                			           N, p_mom, w_par,
		                			           0, 0,
		                			           pot_ptr_nn, pot_ptr_np);
		for (int i = 0; i<1; i++){
			for (int j = 0; j<N; j++){
				//printf("%.5e %.5e %.5e\n", p_par[i], p_par[j], hbarc*MN*V_1S0_interpolated[i*N+j]);
				printf("%.5e %.5e %.5e\n", p_par[i], p_par[j], hbarc*MN*t_unco_array[i*(N+1)+j]);
			}
		}
		return 0;*/

		cout << "Starting Faddeev Iterator" << endl;
		calculate_faddeev_convergence(state_3N_asym_array,
									  P123_array,
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

