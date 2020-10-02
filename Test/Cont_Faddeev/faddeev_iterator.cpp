
#include "faddeev_iterator.h"

double p1(double q, double qp, double x){
    return sqrt(0.25*q*q + qp*qp + x*q*qp); 
}

double p2(double q, double qp, double x){
    return sqrt(q*q + 0.25*qp*qp + x*q*qp);
}

double Atilde (int alpha, int alphaprime, int Ltotal, int Jj_dim, int *L12_Jj, int *l3_Jj, int *J12_Jj, int *two_j3_Jj, int *S12_Jj, int *T12_Jj, int two_J, int two_T, double *SixJ_array, int two_jmax_SixJ){
    double ret = 0;
    double root_term = 0;
    for (int two_Stotal = 1; two_Stotal <= 3; two_Stotal += 2){
        root_term = (2*J12_Jj[alpha] + 1)      * (two_j3_Jj[alpha] + 1)      * (2*S12_Jj[alpha] + 1)      * (2*T12_Jj[alpha] + 1)
                  * (2*J12_Jj[alphaprime] + 1) * (two_j3_Jj[alphaprime] + 1) * (2*S12_Jj[alphaprime] + 1) * (2*T12_Jj[alphaprime] + 1);
                  
        ret +=   sqrt(root_term)
               * gsl_sf_pow_int(-1, S12_Jj[alphaprime] + T12_Jj[alphaprime])
               * (two_Stotal + 1)
               * gsl_sf_coupling_9j(2 * L12_Jj[alpha],     2 * S12_Jj[alpha],     2 * J12_Jj[alpha],     2 * l3_Jj[alpha],     1, two_j3_Jj[alpha],     2 * Ltotal, two_Stotal, two_J)
               * gsl_sf_coupling_9j(2 * L12_Jj[alphaprime], 2 * S12_Jj[alphaprime], 2 * J12_Jj[alphaprime], 2 * l3_Jj[alphaprime], 1, two_j3_Jj[alphaprime], 2 * Ltotal, two_Stotal, two_J)
               * SixJSymbol(SixJ_array, two_jmax_SixJ, 1, 1, 2 * S12_Jj[alpha], 1, two_Stotal, 2 * S12_Jj[alphaprime])
               * SixJSymbol(SixJ_array, two_jmax_SixJ, 1, 1, 2 * T12_Jj[alpha], 1, two_T, 2 * T12_Jj[alphaprime]);

    }
    return ret;
}

double Gtilde_new (double p, double q, double x, int alpha, int alphaprime, int N_alpha, int Lmax, int *L12_Jj, int *l3_Jj, double *A_store, int two_Jtotal){
    double ret = 0.0;
    double fac1, fac2;

    double pi1 = pi1_tilde(p, q, x);
    double pi2 = pi2_tilde(p, q, x);

    double costheta1 = -(0.5 * p + 0.75 * q * x) / pi1;
    double costheta2 = (p - 0.5 * q * x) / pi2;

    int L12 = L12_Jj[alpha];
    int l3 = l3_Jj[alpha];
    int L12prime = L12_Jj[alphaprime];
    int l3prime = l3_Jj[alphaprime];

    for (int Ltotal = max(abs(L12 - l3), abs(L12prime - l3prime)); Ltotal <= min((two_Jtotal + 5) / 2, min(L12 + l3, L12prime + l3prime)); Ltotal++){
        fac1 = 8.0 * M_PI * M_PI * A_store[alpha * N_alpha * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal];

        for (int Mtotal = -min(l3, Ltotal); Mtotal <= min(l3, Ltotal); Mtotal++){
            fac2 = ClebschGordan(2 * L12, 2 * l3, 2 * Ltotal, 0, 2 * Mtotal, 2 * Mtotal)
                   * sqrt((2.0 * L12 + 1) / (4 * M_PI))
                   * gsl_sf_pow_int(-1, Mtotal)
                   * Plm(l3, Mtotal, x); // -1^M phase since azimutal angles of p' and q' = pi

            for (int M12primesum = -L12prime; M12primesum <= L12prime; M12primesum++){
                if (abs(Mtotal - M12primesum) <= l3prime){
                    ret += fac1
                           * fac2
                           * ClebschGordan(2 * L12prime, 2 * l3prime, 2 * Ltotal, 2 * M12primesum, 2 * Mtotal - 2 * M12primesum, 2 * Mtotal)
                           * Plm(L12prime, M12primesum, costheta1)
                           * Plm(l3prime, Mtotal - M12primesum, costheta2);
                }
            }
        }
    }

    return ret;
}

void iterate_faddeev(double** state_3N_symm_array,
                     int& Np, double** p_array, double** wp_array,
                     int& Nq, double** q_array, double** wq_array,
                     int& Nalpha, 
                     int** L_2N_array,
                     int** S_2N_array,
                     int** J_2N_array,
                     int** T_2N_array,
                     int** l_3N_array,
                     int** two_j_3N_array){
    
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
    int l_interpolate_max = l_interpolate_max = 2 * (lmax - 3) + 3;

    cout << "lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";

    int Lmax = max_L12 + max_l3;
    int kLegendremax = 2 * max_L12;

    int two_jmax_Clebsch = 2 * lmax;
    int jmax_Clebsch = lmax;
    int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??

    /* for angular integration in Gtilde */
    double x_Gtilde[Nx_Gtilde];
    double wx_Gtilde[Nx_Gtilde];

    calc_gauss_points (x_Gtilde, wx_Gtilde, -1.0, 1.0, Nx_Gtilde);

    MKL_INT64 g_N = (kLegendremax + 1) * (Lmax + 1) * (Lmax + 1) * Jj_dim * Jj_dim;
    double *gtilde_array = new double[g_N];

    MKL_INT64 Gtilde_N = Jj_dim * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde;
    MKL_INT64 Atilde_N = Jj_dim * Jj_dim * (Lmax + 1);

    double *Gtilde_store = new double[Gtilde_N];
    double *Atilde_store = new double[Atilde_N];

    for (MKL_INT64 i = 0; i <= Atilde_N - 1; i++)
    {
        Atilde_store[i] = 0.0;
    }

    /* Pre-store Wigner-6j symbols */
    double *SixJ_array = new double[(two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)];
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

    pmax_3N = 0.0;
    qmax_3N = 0.0;

    /* Pre-store terms of g that purely depend on alpha and alphaprime */
    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
        for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
            for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
                Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
            }
        }
    }

    /* Pre-store G function as a an array for use in Faddeev equations */
    long int fullsize = Np_3N*Nq_3N*Nx_Gtilde*Jj_dim*(Jj_dim - 1);
    long int counter = 0;
    int frac_n, frac_o=0;
    #pragma omp parallel{
        #pragma omp for
        for (MKL_INT64 p_index = 0; p_index <= Np_3N - 1; p_index++){
            for (MKL_INT64 q_index = 0; q_index <= Nq_3N - 1; q_index++){
                for (MKL_INT64 x_index = 0; x_index <= Nx_Gtilde - 1; x_index++){
                    for (MKL_INT64 alpha = 0; alpha <= Jj_dim - 1; alpha++){
                        for (MKL_INT64 alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
                            Gtilde_store[alpha * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde + alphaprime * Np_3N * Nq_3N * Nx_Gtilde + p_index * Nq_3N * Nx_Gtilde + q_index * Nx_Gtilde + x_index]
                                = Gtilde_new (p_3N[p_index], q_3N[q_index], x_Gtilde[x_index], alpha, alphaprime, Jj_dim, Lmax, L12_Jj, l3_Jj, Atilde_store, two_J);
                            
                            counter += 1;
                            frac_n = (100*counter)/fullsize;
                            if (frac_n>frac_o){cout << frac_n << "%" << endl; frac_o=frac_n;}
                        }
                    }
                }
            }
        }
    }
}