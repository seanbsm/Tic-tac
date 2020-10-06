
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

    /* Free memory */
    delete [] Atilde_store;
    delete [] SixJ_array;
}

void iterate_faddeev(double* state_3N_symm_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                     int idx_alpha_proj, int idx_p_proj, int idx_q_proj,
                     int two_T, int two_J, int PAR){
    
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
    S_spline(S_p_array, p_array, Np);
    S_spline(S_q_array, q_array, Nq);
    /* Spline objects for Faddeev integral interpolations */
    double interpolate_p = 0;
    double interpolate_q = 0;

    /* Construct G_array
     * BEWARE: this is a multidimensional array that can consume a lot of memory
     * (e.g. Nalpha=50, Np=Nq=30, Nx=20: dimension = 4.5*1e7) */
    double* G_array = new double [Nalpha * Nalpha * Np * Nq * Nx];
    make_G_array(G_array, Np, p_array, Nq, q_array, Nx, x_array,
                 Nalpha, L_2N_array, S_2N_array, J_2N_array, T_2N_array, l_3N_array, two_j_3N_array,
                 two_T, two_J);

    /* First faddeev state summation start */

    int    L_p   = 0;
    int    L_pp  = 0;
    int    G_idx = 0;
    double G     = 0;
    double total_sum = 0;
    /* Loop over alpha prime summation */
    for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){
        L_p = L_2N_array[idx_alpha_p];
        /* Loop over alpha double-prime summation */
        for (int idx_alpha_pp=0; idx_alpha_pp<Nalpha; idx_alpha_pp++){
            L_pp = L_2N_array[idx_alpha_pp];
            /* Loop over q prime integral */
            for (int idx_q_p=0; idx_q_p<Nq; idx_q_p++){
                /* Loop over x integral*/
                for (int idx_x=0; idx_x<Nx; idx_x++){
                    
                    /* Evaluate p1 and p2 functions */
                    double pi1 = pi1_tilde(p_array[p_index], q_array[q_index], x_array[idx_x]);
                    double pi2 = pi2_tilde(p_array[p_index], q_array[q_index], x_array[idx_x]);
                    
                    /* Determine which spline-array index we need */
                    int pi1_index = 0; while ((p_array[pi1_index] < pi1) && (pi1_index < Np - 1)) pi1_index++; if (pi1_index > 0) pi1_index--;
                    int pi2_index = 0; while ((q_array[pi2_index] < pi2) && (pi2_index < Nq - 1)) pi2_index++; if (pi2_index > 0) pi2_index--;

                    /* Calculate spline functions from pre-calculated spline-arrays */
                    /* If-test to avoid artefacts at the upper limits of mesh system */
                    if (pi1 <= p_array[Np - 1]){
                        interpolate_p = (
                                            S_p_array[pi1_index * Np + pprime_index]
                                            + (pi1 - p_array[pi1_index]) * S_p_array[pi1_index * Np + pprime_index + Np * (Np - 1)]
                                            + (pi1 - p_array[pi1_index]) * (pi1 - p_array[pi1_index]) * S_p_array[pi1_index * Np + pprime_index + 2 * Np * (Np - 1)]
                                            + (pi1 - p_array[pi1_index]) * (pi1 - p_array[pi1_index]) * (pi1 - p_array[pi1_index]) * S_p_array[pi1_index * Np + pprime_index + 3 * Np * (Np - 1)]
                                        );
                    }
                    else{
                        interpolate_p = 0.0;
                    }

                    /* If-test to avoid spline artefacts at the upper limits of mesh system */
                    if (pi2 <= q_array[Nq - 1]){
                        interpolate_q = (
                                            S_q_array[pi2_index * Nq + qprime_index]
                                            + (pi2 - q_array[pi2_index]) * S_q_array[pi2_index * Nq + qprime_index + Nq * (Nq - 1)]
                                            + (pi2 - q_array[pi2_index]) * (pi2 - q_array[pi2_index]) * S_q_array[pi2_index * Nq + qprime_index + 2 * Nq * (Nq - 1)]
                                            + (pi2 - q_array[pi2_index]) * (pi2 - q_array[pi2_index]) * (pi2 - q_array[pi2_index]) * S_q_array[pi2_index * Nq + qprime_index + 3 * Nq * (Nq - 1)]
                                        );
                    }
                    else{
                        interpolate_q = 0.0;
                    }
                    
                    /* Retrieve G_{alpha,alpha''}(q,q'',x) from pre-calculated G_array */
                    G_idx = idx_alpha_proj * Nalpha*Np*Nq*Nx + idx_alpha_pp * Np*Nq*Nx + idx_q_proj * Nq*Nx + idx_q_p * Nx + idx_x;
                    G = G_array[G_idx];

                    /* Evaluate integral term and add to total sum */
                    total_sum += wq_array[idx_q_p]*q_array[idx_q_p]*q_array[idx_q_p] * wx_array[idx_x] * t_element * G * interpolate_q * prev_psi;
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