#include "make_swp_states.h"

/* Finds the eigenvalues and eigenvectors
 * of a real, symmetric matrix A.
 * For simplicity here, A must be
 * stored as an upper triangle. w will
 * be filled with the eigenvalues in
 * ascending order, and z will be a matrix
 * with the corresponding eigenvectors
 * (column by column). Lastly, we work
 * with row major matrices, as is usual
 * with C & C++ */
 
void diagonalize_real_symm_matrix(float *A, float *w, float *z, int N){
	char jobz = 'V';
	char uplo = 'U';
	
	LAPACKE_sspevd(LAPACK_ROW_MAJOR, jobz, uplo, N, A, w, z, N);
}

void diagonalize_real_symm_matrix(double *A, double *w, double *z, int N){
	char jobz = 'V';
	char uplo = 'U';
	
	LAPACKE_dspevd(LAPACK_ROW_MAJOR, jobz, uplo, N, A, w, z, N);
}

void construct_free_hamiltonian(double* H0_WP_array,
							    int Np_WP, double* p_WP_array){
	
    /* Loop over p-momenta */
    for (int idx_p=0; idx_p<Np_WP; idx_p++){
		/* Kinetic energy at bin-boundaries */
		double p1 = p_WP_array[idx_p];
		double p2 = p_WP_array[idx_p+1];
	
		/* Free Hamiltonian for 2-nucleon pair for momentum WPs */
		H0_WP_array[idx_p] = (p2*p2 + p2*p1 + p1*p1)/(3*MN);
		/* Free Hamiltonian for 2-nucleon pair for energy WPs */
        //H0_WP_array[idx_p] = (p2*p2 + p1*p1)/(2*MN);
	}
}

/* Constructs a NN-pair full Hamiltonians as upper-triangular from given arrays */
void construct_full_hamiltonian(double* mat_ptr_H,
								double* mat_ptr_H0,
								double* mat_ptr_V,
                                int     mat_dim){

    /* Row p-momentum index loop */
    for (int idx_bin_r=0; idx_bin_r<mat_dim; idx_bin_r++){
					
		/* Because we only store the upper triangular part we need special indexing for H */
		int idx_r_r = mat_dim*(mat_dim-1)/2 - (mat_dim-idx_bin_r)*(mat_dim-idx_bin_r-1)/2 + idx_bin_r;

		/* Set diagonal to free kinetic energy given by H0_array */
		mat_ptr_H[idx_r_r] += mat_ptr_H0[idx_bin_r];

        /* Column p-momentum index loop (note it starts on idx_bin_r since H is upper triangular) */
        for (int idx_bin_c=idx_bin_r; idx_bin_c<mat_dim; idx_bin_c++){
						
			/* Because we only store the upper triangular part we need special indexing for H */
			int idx_r_c = mat_dim*(mat_dim-1)/2 - (mat_dim-idx_bin_r)*(mat_dim-idx_bin_r-1)/2 + idx_bin_c;

			/* Add potential element to Hamiltonian */
			mat_ptr_H[idx_r_c] += mat_ptr_V[idx_bin_r*mat_dim + idx_bin_c];
		}
	}
}

/* This function reorders the eigenvalues and corresponding vectors to correspond to coupled channels.
 * To understand why it's best to read up on WPCD theory - it's too long to explain here. */
void reorder_coupled_eigenspectrum(double* eigenvalues,
								   double* eigenvectors,
								   int     Np_WP){
	
	/* Temporary arrays to hold values as we reorder indices */
	double buffer_array_vals [2*Np_WP];
	double buffer_array_vecs [4*Np_WP*Np_WP];

	for (int i=0; i<Np_WP; i++){
		/* Temporarily store reorder of eigenvalues */
		buffer_array_vals[i] 	   = eigenvalues[2*i];
		buffer_array_vals[i+Np_WP] = eigenvalues[2*i+1];

		for (int j=0; j<2*Np_WP; j++){
			/* Temporarily store reorder of eigenvectors */
			buffer_array_vecs[j*2*Np_WP + i] 	     = eigenvectors[j*2*Np_WP + 2*i];
			buffer_array_vecs[j*2*Np_WP + i + Np_WP] = eigenvectors[j*2*Np_WP + 2*i+1];
		}
	}

	/* Write reordered elements back into original arrays */
	for (int i=0; i<2*Np_WP; i++){
		eigenvalues[i] = buffer_array_vals[i];
		for (int j=0; j<2*Np_WP; j++){
			eigenvectors[i*2*Np_WP + j] = buffer_array_vecs[i*2*Np_WP + j];
		}
	}
}

void look_for_unphysical_bound_states(double* eigenvalues,
						   			  int     mat_dim,
						   			  bool    chn_3S1){
											 
	int num_bound_states_found = 0;

	/* Cound number of bound states in eigenspectrum */
	for (int idx=0; idx<mat_dim; idx++){
		if (eigenvalues[idx]<0){
			num_bound_states_found += 1;
		}
		//if (chn_3S1){
		//	std::cout << eigenvalues[idx] << std::endl;
		//}
	}
	
	/* See if we find the expected number of bounds states (1 for 3S1, 0 otherwise),
	 * if not then we abort program (unphysical scenario) */
	
	if (num_bound_states_found==0 and chn_3S1==false){}
	else if (num_bound_states_found==1 and chn_3S1==true){}
	else{
		raise_error("Found an unphysical bound state in NN-pair Hamiltonian eigenspectrum!");
	}
}

void make_swp_bin_boundaries(double* eigenvalues,
							 double* e_SWP_array,
							 int	 Np_WP,
							 bool    coupled,
							 bool    chn_3S1){

	/* Set first boundary */
	e_SWP_array[0] = 0;
	if (coupled){
		e_SWP_array[Np_WP] = 0;
	}
	
	/* Set all the middle boundaries */
	for (int idx_p=0; idx_p<Np_WP-1; idx_p++){
		e_SWP_array[idx_p+1] = 0.5*(eigenvalues[idx_p+1] + eigenvalues[idx_p]);

		if (coupled){
			e_SWP_array[idx_p+1+Np_WP] = 0.5*(eigenvalues[idx_p+1+Np_WP] + eigenvalues[idx_p+Np_WP]);
		}
	}

	/* Set end boundary */
	e_SWP_array[Np_WP] = eigenvalues[Np_WP-1] + 0.5*(eigenvalues[Np_WP-2] + eigenvalues[Np_WP-3]);
	if (coupled){
		e_SWP_array[2*Np_WP] = eigenvalues[2*Np_WP-1] + 0.5*(eigenvalues[2*Np_WP-2] + eigenvalues[2*Np_WP-3]);
	}

	/* Modify first boundary for 3S1 */
	if (chn_3S1){
		e_SWP_array[0] = eigenvalues[0];
		e_SWP_array[1] = 0;
	}
}

void make_swp_states(double* e_SWP_unco_array,
					 double* e_SWP_coup_array,
					 double* C_WP_unco_array,
					 double* C_WP_coup_array,
					 double* V_WP_unco_array,
                     double* V_WP_coup_array,
					 int Np_WP, double* p_WP_array,
					 int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array,
					 int J_2N_max){
	
	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = 2*(J_2N_max+1);
	int num_coup_chns =    J_2N_max;

	/* Check-lists to keep track of which 2N Hamiltonian diagonalizations have been done.
	 * This removes excessive work due to distinct 3N channels containing equal
	 * 2N channels, as well as overwriting existing calculations (thus giving wrongful results) */
	bool* check_list_unco = new bool [num_unco_chns];
	bool* check_list_coup = new bool [num_coup_chns];

	/* Set check_list-arrays to false */
	for (int i=0; i<num_unco_chns; i++){
		check_list_unco[i] = false;
	}
	for (int i=0; i<num_coup_chns; i++){
		check_list_coup[i] = false;
	}

	/* Boundaries of scattering wave-packets (SWP) in energy */
	int e_SWP_unco_array_size =   (Np_WP+1) * num_unco_chns;
	int e_SWP_coup_array_size = 2*(Np_WP+1) * num_coup_chns;

	/* Set e_SWP-arrays to zero */
	for (int i=0; i<e_SWP_unco_array_size; i++){
		e_SWP_unco_array[i] = 0;
	}
	for (int i=0; i<e_SWP_coup_array_size; i++){
		e_SWP_coup_array[i] = 0;
	}

	/* The Hamiltonian matrices are symmetric and stored as upper-triangular
	 * for computational efficiency when diagonalizing */
	int H_unco_array_size =   Np_WP*(  Np_WP+1)/2 * num_unco_chns;
    int H_coup_array_size = 2*Np_WP*(2*Np_WP+1)/2 * num_coup_chns;
	double* H_WP_unco_array = new double [H_unco_array_size];
	double* H_WP_coup_array = new double [H_coup_array_size];

	/* Set H_WP-arrays to zero */
	for (int i=0; i<H_unco_array_size; i++){
		H_WP_unco_array[i] = 0;
	}
	for (int i=0; i<H_coup_array_size; i++){
		H_WP_coup_array[i] = 0;
	}

	/* Free Hamiltonian arrays */
	int H0_unco_array_size =   Np_WP;
    int H0_coup_array_size = 2*Np_WP;
	double* H0_WP_unco_array = new double [H0_unco_array_size];
	double* H0_WP_coup_array = new double [H0_coup_array_size];

	/* Set H0_WP-arrays to zero */
	for (int i=0; i<H0_unco_array_size; i++){
		H0_WP_unco_array[i] = 0;
	}
	for (int i=0; i<H0_coup_array_size; i++){
		H0_WP_coup_array[i] = 0;
	}

	/* Fill free Hamiltonian arrays. All 3 calls will produce the same,
	 * but it's quite fast and makes the code clearer to interpret */
	construct_free_hamiltonian( H0_WP_unco_array,        Np_WP, p_WP_array);
	construct_free_hamiltonian(&H0_WP_coup_array[0],     Np_WP, p_WP_array);
	construct_free_hamiltonian(&H0_WP_coup_array[Np_WP], Np_WP, p_WP_array);
	
	/* Row state */
    for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        int L_r = L_2N_array[idx_alpha_r];
        int S_r = S_2N_array[idx_alpha_r];
        int J_r = J_2N_array[idx_alpha_r];
        int T_r = T_2N_array[idx_alpha_r];
        
        /* Column state */
        for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
            int L_c = L_2N_array[idx_alpha_c];
            int S_c = S_2N_array[idx_alpha_c];
            int J_c = J_2N_array[idx_alpha_c];
            int T_c = T_2N_array[idx_alpha_c];

            /* Check if possible channel through interaction */
            if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){

                /* Detemine if this is a coupled channel */
	            bool coupled_matrix = false;
	            if (L_r!=L_c or (L_r==L_c and L_r!=J_r and J_r!=0)){ // This counts 3P0 as uncoupled (which is intended)
	            	coupled_matrix  = true;
	            }

                /* Skip redundant calculations by only doing the coupled calculation when L_r<L_c */
                if (coupled_matrix){
                    if ( (L_r<L_c)==false ){
                        continue;
                    }
                }

				/* Hamiltonian matrix pointer and dimension
                 * Indexing format of Hamitonian arrays: (channel index)*(num rows)*(num columns) + (row index)*(row length) + (column index) */
				int     mat_dim   		= 0;
				int		chn_idx	  		= 0;
				double* mat_ptr_H 		= NULL;
				double* mat_ptr_H0 		= NULL;
				double* mat_ptr_V 		= NULL;
				double* mat_ptr_C 		= NULL;
				double* e_SWP_array_ptr = NULL;
                if (coupled_matrix){
					mat_dim   		= 2*Np_WP;
					chn_idx   		= J_r-1;

					/* Check if 2N channels diagonalization has already
					 * been performed in previous loop-iterations,
					 * and if not then set to true */
					if (check_list_coup[chn_idx]==true){
						continue;
					}
					else{
						check_list_coup[chn_idx] = true;
					}

					/* H-matrices are stored as upper-triangular -> special indexing and step-length */
					mat_ptr_H 		= &H_WP_coup_array [chn_idx * mat_dim*(mat_dim+1)/2];
					mat_ptr_V 		= &V_WP_coup_array [chn_idx * mat_dim*mat_dim];
					mat_ptr_C 		= &C_WP_coup_array [chn_idx * mat_dim*mat_dim];
					e_SWP_array_ptr = &e_SWP_coup_array[chn_idx * mat_dim];
					mat_ptr_H0 		= H0_WP_coup_array;
                }
			    else{
					mat_dim   		= Np_WP;
					chn_idx   		= 2*J_r + S_r;

					/* Check if 2N channels diagonalization has already
					 * been performed in previous loop-iterations,
					 * and if not then set to true */
					if (check_list_unco[chn_idx]==true){
						continue;
					}
					else{
						check_list_unco[chn_idx] = true;
					}

					/* H-matrices are stored as upper-triangular -> special indexing and step-length */
					mat_ptr_H 		= &H_WP_unco_array [chn_idx * mat_dim*(mat_dim+1)/2];
					mat_ptr_V 		= &V_WP_unco_array [chn_idx * mat_dim*mat_dim];
					mat_ptr_C 		= &C_WP_unco_array [chn_idx * mat_dim*mat_dim];
					e_SWP_array_ptr = &e_SWP_unco_array[chn_idx * mat_dim];
					mat_ptr_H0 		= H0_WP_unco_array;
                }

				/* Construct channel Hamiltonian */
				construct_full_hamiltonian(mat_ptr_H,
										   mat_ptr_H0,
										   mat_ptr_V,
                                		   mat_dim);

				/* Hamiltonian eigenvalue array */
				double eigenvalues [mat_dim];

				/* Diagonalize channel Hamiltonian - fill eigenvalues and C_array coefficients */
				diagonalize_real_symm_matrix(mat_ptr_H, eigenvalues, mat_ptr_C, mat_dim);

				/* Abort if unphysical bound states are found in eigenvalues,
				   or if 3S1-bound state is missing */
				bool chn_3S1 = (coupled_matrix && J_r==1);
				look_for_unphysical_bound_states(eigenvalues, mat_dim, chn_3S1);
				
				/* The eigenspectrum of coupled channels is returned in ascending
				 * order from the diagonalization routine, so we reorder the spectrum here */
				if (coupled_matrix){
					reorder_coupled_eigenspectrum(eigenvalues,
												  mat_ptr_C,
												  Np_WP);
				}
				
				/* Construct energy bin boundaries for swp states */
				make_swp_bin_boundaries(eigenvalues,
										e_SWP_array_ptr,
										Np_WP,
										coupled_matrix,
										chn_3S1);
			}
		}
	}

	/* Delete temporary arrays */
	delete [] check_list_unco;
	delete [] check_list_coup;

	delete [] H_WP_unco_array;
	delete [] H_WP_coup_array;

	delete [] H0_WP_unco_array;
	delete [] H0_WP_coup_array;
}