
#include "faddeev_iterator.h"

/* Copied from https://stackoverflow.com/questions/25313816/gram-schmidt-function-not-working-c/25320243
 * The function uses the modified Gram-Scmidt orgonalization routine to read a matrix of vectors (state_matrix)
 * to create an orgothogonal basis from them (state_basis). */
void modified_gram_schmidt(double* state_matrix, double* state_basis, int N){
    double r [N*N];
    double v [N*N];

    for (int i=0; i < N; i++){
        for (int j=0; j<N; j++){
            //v[i][j] = matrix[i][j];
            v[i*N + j] = state_matrix[i*N + j];
        }
    }

    for (int i=0; i<N; i++){
        //r[i][i] = getNorm(v[i]);
        r[i*N + i] = 0;
        for (int n=0; n<N; n++){
            r[i*N + i] += v[i*N + n] * v[i*N + n];
        }
        r[i*N + i] = sqrt(r[i*N + i]);

        for (int j=0; j<N; j++){
            //base[i][j] = v[i][j] / r[i][i];
            state_basis[i*N + j] = v[i*N + j] / r[i*N + i];
        }

        for (int k=i+1;  k<N; k++){
            //r[i][k] = dotProduct(base[i],v[k]);
            r[i*N + k] = 0;
            for (int n=0; n<N; n++){
                r[i*N + k] += state_basis[i*N + n] * v[k*N + n];
            }
            
            for (int j=0; j<N; j++){
                //v[k][j] = v[k][j] - r[i][k] * base[i][j];
                v[k*N + j] = v[k*N + j] - r[i*N + k] * state_basis[i*N + j];
            }
        }
    }
}

void make_G_array(double* G_array,
                  int Np, double* p_array,
                  int Nq, double* q_array,
                  int Nx, double* x_array,
                  int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                  int two_T, int two_J){
    
    /* Change in notation from my program to script in this function */
    int     Jj_dim       = Nalpha;
    int     Nx_Gtilde    = Nx;
    int     Np_3N        = Np;
    int     Nq_3N        = Nq;
    int*    L12_Jj       = l_3N_array;
    int*    l3_Jj        = L_2N_array;
    int*    J12_Jj       = J_2N_array;
    int*    two_j3_Jj    = two_j_3N_array;
    int*    S12_Jj       = S_2N_array;
    int*    T12_Jj       = T_2N_array;
    double* p_3N         = p_array;
    double* q_3N         = q_array;
    double* x_Gtilde     = x_array;
    double* Gtilde_store = G_array;
    /* The remaining function is copied from an external script and utilises auxiliary.cpp functionality
     * (I have commented out print-commands) */
    
    // determine optimized Lmax: Lmax = max(get_L)+max(get_l)
    int max_L12 = 0;
    int max_l3 = 0;
    int max_J12 = 0;

    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++)
    {
        if (J12_Jj[alpha] > max_J12) max_J12 = J12_Jj[alpha];
        if (L12_Jj[alpha] > max_L12) max_L12 = L12_Jj[alpha];
        if (l3_Jj[alpha] > max_l3) max_l3 = l3_Jj[alpha];
    }

    // for F_local_matrix prestorage and F_interpolate
    int lmax = GSL_MAX_INT(max_l3, max_L12) + 3; // for C4 it is possible to couple l=lmax with THREE Y_{1}^{mu}
    //int l_interpolate_max = 2 * (lmax - 3) + 3;

    //cout << "lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";

    int Lmax = max_L12 + max_l3;
    //int kLegendremax = 2 * max_L12;

    //int two_jmax_Clebsch = 2 * lmax;
    //int jmax_Clebsch = lmax;
    int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??


    // for angular integration in Gtilde
    //double x_Gtilde[Nx_Gtilde];
    //double wx_Gtilde[Nx_Gtilde];

    //calc_gauss_points (x_Gtilde, wx_Gtilde, -1.0, 1.0, Nx_Gtilde);

    //MKL_INT64 g_N = (kLegendremax + 1) * (Lmax + 1) * (Lmax + 1) * Jj_dim * Jj_dim;
    //double *gtilde_array = new double[g_N];

    //MKL_INT64 Gtilde_N = Jj_dim * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde;
    MKL_INT64 Atilde_N = Jj_dim * Jj_dim * (Lmax + 1);

    //double *Gtilde_store = new double[Gtilde_N];
    double *Atilde_store = new double[Atilde_N];

    for (MKL_INT64 i = 0; i <= Atilde_N - 1; i++){
        Atilde_store[i] = 0.0;
    }

    //cout << "prestore SixJ...\n";
    double *SixJ_array = new double[(two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)];
    //cout << "SixJ_test (>0?):" << (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) << "\n";

    std::cout << "   - Pre-storing Wigner-6j coupling symbols ..." << std::endl;
    #pragma omp parallel for collapse(3)
    
    for (int two_l1 = 0; two_l1 <= two_jmax_SixJ; two_l1++){
        for (int two_l2 = 0; two_l2 <= two_jmax_SixJ; two_l2++){
            for (int two_l3 = 0; two_l3 <= two_jmax_SixJ; two_l3++){
                for (int two_l4 = 0; two_l4 <= two_jmax_SixJ; two_l4++){
                    for (int two_l5 = 0; two_l5 <= two_jmax_SixJ; two_l5++){
                        for (int two_l6 = 0; two_l6 <= two_jmax_SixJ; two_l6++){
                            SixJ_array[
                                two_l1 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l2 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l3 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l4 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
                                + two_l5 * (two_jmax_SixJ + 1)
                                + two_l6
                            ]
                            // implement checks for quantum numbers because of bug in gsl library
                                = gsl_sf_coupling_6j(two_l1, two_l2, two_l3, two_l4, two_l5, two_l6);
                        }
                    }
                }
            }
        }
    }

    std::cout << "   - Calculating Atilde ..." << std::endl;
    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
        for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
            for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
                Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
            }
        }
    }

    /* Calculate G_array */
    //long int fullsize = Np_3N*Nq_3N*Nx_Gtilde*Jj_dim*(Jj_dim - 1);
    //long int counter = 0;
    //int frac_n, frac_o=0;
    std::cout << "   - Calculating Gtilde ..." << std::endl;
    #pragma omp parallel
    {
        #pragma omp for

        for (MKL_INT64 p_index = 0; p_index <= Np_3N - 1; p_index++){
            for (MKL_INT64 q_index = 0; q_index <= Nq_3N - 1; q_index++){
                for (MKL_INT64 x_index = 0; x_index <= Nx_Gtilde - 1; x_index++){
                    for (MKL_INT64 alpha = 0; alpha <= Jj_dim - 1; alpha++){
                        for (MKL_INT64 alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
                            Gtilde_store[alpha * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde + alphaprime * Np_3N * Nq_3N * Nx_Gtilde + p_index * Nq_3N * Nx_Gtilde + q_index * Nx_Gtilde + x_index]
                                = Gtilde_new (p_3N[p_index], q_3N[q_index], x_Gtilde[x_index], alpha, alphaprime, Jj_dim, Lmax, L12_Jj, l3_Jj, Atilde_store, two_J);
                            
                            //counter += 1;
                            //frac_n = (100*counter)/fullsize;
                            //if (frac_n>frac_o){cout << frac_n << "%" << endl; frac_o=frac_n;}
                        }
                    }
                }
            }
        }
    }
    std::cout << "   - Done!" << std::endl;
    /* Free memory */
    delete [] Atilde_store;
    delete [] SixJ_array;
}

/* Indexing: S[j*Nx + i + k*Nx(Nx-1)] for S_i(x) where x is in bin j, and k is the coefficient we need
 * We have a function evaluation at x_i. 
 * This means we use x_j when x is in bin j, i.e. (x-x_j). */
double interpolate_using_spline_array(double* S_array, int Np, double* p_array, int idx_pi, double p, int idx_p){
    double interpolate = 0;

    double pi = p_array[idx_pi];

    if (p <= p_array[Np - 1]){
        interpolate = (
                                                           S_array[idx_pi * Np + idx_p]
                            + (p - pi)                   * S_array[idx_pi * Np + idx_p +     Np * (Np - 1)]
                            + (p - pi)*(p - pi)          * S_array[idx_pi * Np + idx_p + 2 * Np * (Np - 1)]
                            + (p - pi)*(p - pi)*(p - pi) * S_array[idx_pi * Np + idx_p + 3 * Np * (Np - 1)]
                        );
    }

    return interpolate;
}

void iterate_faddeev(double* state_3N_symm_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                     int idx_alpha, int idx_p, int idx_q,
                     int two_T, int two_J, int PAR,
                     potential_model* pot_ptr_np,
                     potential_model* pot_ptr_nn){
    
    int Faddeev_iterations = 10;

    /* Matrix of states to be orthogonalized using Gram-Schmidt algorithm
     * and the nusing Lanczos algorithm to solve Faddeev */
    double* psi_matrix = new double [Nalpha*Np*Nq * Faddeev_iterations];

    /* Make x integral mesh */
    int Nx = 20;
    double* x_array  = new double [Nx];
    double* wx_array = new double [Nx];
    calc_gauss_points (x_array, wx_array, -1., 1., Nx);

    /* Declare spline arrays */
    double* S_p_array = new double [4*Np*(Np-1)];
    double* S_q_array = new double [4*Nq*(Nq-1)];
    /* Calculate and fill spline arrays
     * In steps of 4, these arrays contain the coefficients for quadratic interpolation of momenta.
     * See Eq. (II.16) in https://doi.org/10.1007/BF01417437 */
     std::cout << "Calculating spline-objects" << std::endl;
    S_spline(S_p_array, p_array, Np);
    S_spline(S_q_array, q_array, Nq);
    std::cout << "   - Done!" << std::endl;

    /* Construct G_array
     * BEWARE: this is a multidimensional array that can consume a lot of memory
     * (e.g. Nalpha=50, Np=Nq=30, Nx=20: dimension = 4.5*1e7) */
    std::cout << "Calculating G_array" << std::endl;
    double* G_array = new double [Nalpha * Nalpha * Np * Nq * Nx];
    make_G_array(G_array, Np, p_array, Nq, q_array, Nx, x_array,
                 Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array, l_3N_array, two_j_3N_array,
                 two_T, two_J);

    /* Current Faddeev iteration */
    int    current_iteration = 0;

    /* Projection state quantum numbers (alpha) */
    int L     = L_2N_array[idx_alpha];
    int S     = S_2N_array[idx_alpha];
    int J     = J_2N_array[idx_alpha];
    int T     = T_2N_array[idx_alpha];
    int two_l = l_3N_array[idx_alpha];
    int two_j = two_j_3N_array[idx_alpha];

    /* Triton ground-state energy (to be determined through Lanczos iterations) */
    double E = 0;
    /* Upper-case is for the 2N-pair, lower-case is for orbiting spectator */
    int L_p  = 0, S_p  = 0, J_p  = 0, T_p  = 0, two_l_p  = 0, two_j_p  = 0;
    int L_pp = 0, S_pp = 0, J_pp = 0, T_pp = 0, two_l_pp = 0, two_j_pp = 0;
    /* Spline object for Faddeev integral interpolations */
    double S_interpolation = 0;
    /* Geometric function indexing and function value */
    int    idx_G     = 0;
    double G_element = 0;
    /* T-matrix element of interest */
    double t_element = 0;
    /* Previous value of wave-function, kept and updated in psi_matrix */
    double prev_psi = 0;
    /* Faddeev summation term */
    double total_sum = 0;
    /* Loop over alpha prime summation */
    std::cout << "Starting Faddeev iterations" << std::endl;
    for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){
        std::cout << idx_alpha_p << " / " << Nalpha << std::endl;
        L_p     = L_2N_array[idx_alpha_p];
        S_p     = S_2N_array[idx_alpha_p];
        J_p     = J_2N_array[idx_alpha_p];
        T_p     = T_2N_array[idx_alpha_p];
        two_l_p = l_3N_array[idx_alpha_p];
        two_j_p = two_j_3N_array[idx_alpha_p];
        
        /* Loop over alpha double-prime summation */
        for (int idx_alpha_pp=0; idx_alpha_pp<Nalpha; idx_alpha_pp++){
            L_pp     = L_2N_array[idx_alpha_pp];
            S_pp     = S_2N_array[idx_alpha_pp];
            J_pp     = J_2N_array[idx_alpha_pp];
            T_pp     = T_2N_array[idx_alpha_pp];
            two_l_pp = l_3N_array[idx_alpha_pp];
            two_j_pp = two_j_3N_array[idx_alpha_pp];

            /* The two-body force can only allow for L_p!=L_pp, everything else must be the same */
            if (two_l_p==two_l_pp and two_j_p==two_j_pp and T_p==T_pp and J_p==J_pp and S_p==S_pp and abs(L_p-L_pp)<=2){
                /* Loop over q prime integral */
                for (int idx_q_p=0; idx_q_p<Nq; idx_q_p++){
                    /* Loop over x integral*/
                    for (int idx_x=0; idx_x<Nx; idx_x++){

                        /* Evaluate p1 and p2 functions */
                        double p1 = pi1_tilde(q_array[idx_q], q_array[idx_q_p], x_array[idx_x]);
                        double p2 = pi2_tilde(q_array[idx_q], q_array[idx_q_p], x_array[idx_x]);
                        //std::cout << "1" << std::endl;

                        /* Determine which spline-array index we need */
                        int idx_p2 = 0; while ((p_array[idx_p2] < p2) && (idx_p2 < Nq - 1)) idx_p2++; if (idx_p2 > 0) idx_p2--;

                        /* Calculate spline functions from pre-calculated spline-arrays */
                        S_interpolation = interpolate_using_spline_array(S_p_array, Np, p_array, idx_p2, p2, idx_p);
                         //std::cout << "2" << std::endl;

                        /* Solve Lippmann-Schwinger equation for current q', q, alpha', and alpha */
                        /* !!! NOTE I USE NUCLEON MASS HERE, THIS IS SLIGHTLY WRONG BUT SHOULD MATTER LITTLE !!! */
                        double E_LS = E - q_array[idx_q_p]*q_array[idx_q_p] / (2 * 3*MN/4);
                        t_element = calculate_t_element(L, L_p, S, J, T,
						                                E_LS, p_array[idx_p], MN,
						                                Np, p_array, wp_array,
						                                idx_p, Nq,
						                                pot_ptr_nn, pot_ptr_np);
                         //std::cout << "3" << std::endl;

                        /* Copy previous Faddeev solution (wave-function) */
                        if (current_iteration==0){  // For the first Faddeev-iteration we use the first column of the matrix
                            prev_psi = psi_matrix[idx_alpha*Nq*Np* + idx_q*Np + idx_p];
                        }
                        else{   // For current_iteration>0 Faddeev-iteration we use the previous column of the matrix
                            prev_psi = psi_matrix[(current_iteration-1)*Nalpha*Nq*Np + idx_alpha*Nq*Np* + idx_q*Np + idx_p];
                        }

                        /* Retrieve G_{alpha,alpha''}(q,q'',x) from pre-calculated G_array */
                        idx_G = idx_alpha * Nalpha*Np*Nq*Nx + idx_alpha_pp * Np*Nq*Nx + idx_q * Nq*Nx + idx_q_p * Nx + idx_x;
                        G_element = G_array[idx_G];
                         //std::cout << "4" << std::endl;

                        double denominator = 1./(pow(p1, L_pp) * pow(p2, L_p));

                        /* Evaluate integral term and add to total sum */
                        total_sum +=  wq_array[idx_q_p] * q_array[idx_q_p] * q_array[idx_q_p]
                                    * wx_array[idx_x] * t_element * G_element * denominator * S_interpolation * prev_psi;
                         //std::cout << "5" << std::endl;
                    }
                }
            }
        }
    }

    /* First faddeev state summation end */


    /* Free memory */
    delete [] x_array;
    delete [] wx_array;
    delete [] S_p_array;
    delete [] S_q_array;
    delete [] G_array;
}