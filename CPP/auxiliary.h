#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm> // std::max, std::min
#include "mkl.h"
#include <fstream>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_fit.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

#include "make_wp_states.h"

using namespace std;

double P (double x, void *params);

double dP (double x, void *params);

void fdf (double x, void *params, double *y, double *dy);

double Legendre (double x, int n);

double dLegendre (double x, int n);

void roots (double *x, int n);

void calc_gauss_points (double *q_par, double *w_par, double pmin_par, double pmax_par, int N_par);

// transform \int_a^infinity f(x) dx = \int_{arctan(a)}^(Pi/2) f(tan(t)) (1+tan(t)^2) dt
// i.e. transform
void calc_gauss_points_tan (double *q_par, double *w_par, double a, double b, int N_par);

// enforce fewer points around the critical areas around x=1 and x=-1
void calc_gauss_points_ang_opt (double *q_par, double *w_par, double pmin_par, double pmax_par, int N_par);

double ClebschGordan (int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);

double ClebschGordan_array (int two_jmax, double *data_array, int two_j1, int two_j2, int two_j3, int two_m1, int two_m2);

double ThreeJarray (int jmax, double *data_array, int j1, int j2, int j3, int m1, int m2);

double SixJSymbol (double *SixJ_array, int two_jmax, int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6);


double NineJSymbol (int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6, int two_j7, int two_j8);

double ClebschGordan_array (int two_jmax, double *data_array, int two_j1, int two_j2, int two_j3, int two_m1, int two_m2);

int KroneckerDelta(int i, int j);

// actually not needed....
complex<double> Ylm (int l, int m, double theta, double phi);

// associate Legendre Polynomials with proper normalization, so Ylm = e^{i m phi} Plm
double Plm (int l, int m, double costheta);

// calc interpolating function: p_par(p_par_index) < x < p_par(p_par_index+1)
// f(x) = \sum_i f(p_par(i)) S_fct(i,x,p_par,N_par,S_vec)
double S_fct(int i, double x, double *p_par, int N_par, double *S_vec_par);

// more efficient method, assumes that j_index is already known, saves the internal loop
double S_fct(int i, double x, int j_index, double *p_par, int N_par, double *S_vec_par);

// cubic spline
void S_spline (double *S, double *p_par, int N_par);

double lin_poly_interpolation (double x1, double y1, double x2, double y2, double x);

complex<double> clin_poly_interpolation (double x1, complex<double> y1, double x2, complex<double> y2, double x);

double quadr_poly_interpolation (double x1, double y1, double x2, double y2, double x3, double y3, double x);

// f1, f2, f3, f4 function values on a rectangle counted counterclockwise beginning in the lower left corner
double lin_poly_2dim (double x1, double x2, double y1, double y2, double f1, double f2, double f3, double f4, double x, double y, double cutoff);

// f1, f2, f3, f4, f5, f6, f7, f8 function values on a rectangle counted counterclockwise beginning in the lower left corner, first rectangle f1-f4 lower coordinate in z, f5-f8 higher coordinate
double lin_poly_3dim (double x1, double x2, double y1, double y2, double z1, double z2,
                      double f1, double f2, double f3, double f4, double f5, double f6, double f7, double f8,
                      double x, double y, double z, double cutoff);

double Floc_l (int l1, int l2, int l1prime, int l2prime,
               int l, double p, double q, double pprime,
               double qprime, int Npq, double *ClebschGordan_data,
               int two_jmax, double *p_tilde_vec, double *wp_tilde_vec,
               double *q_tilde_vec, double *wq_tilde_vec, double *F_interpolate,
               int lmax, double *Plm_p_1, double *Plm_p_2, double *Plm_q_1, double *Plm_q_2);

double Gtilde (double p, double q, double x, int alpha, int alphaprime,
               int Jj_dim, int Lmax, int kLegendremax, double *g_array_par,
               int *L12_Jj, int *l3_Jj);

double pi1_tilde(double p, double q, double x);

double pi2_tilde(double p, double q, double x);

double pi1_prime_tilde(double p, double q, double x);


double pi2_prime_tilde(double p, double q, double x);

double calculate_P123_element_in_WP_basis_mod (double* Gtilde_subarray,
											   int  p_idx_WP, int  q_idx_WP,
											   int pp_idx_WP, int qp_idx_WP,
											   int Np_WP,     double *p_array_WP_bounds,
											   int Nq_WP,     double *q_array_WP_bounds,
											   int Nx, double *x_array, double *wx_array,
											   int Nphi,
											   double* sin_phi_subarray,
											   double* cos_phi_subarray,
											   double* wphi_subarray);

void generate_Ptilde_new(double *P123_store, MKL_INT64 Pdim,
                         MKL_INT64 N_p, double *p,
                         MKL_INT64 N_q, double *q,
                         MKL_INT64 N_x, double *x, double *wx,
                         MKL_INT64 Jj_dim, double pmax, double qmax,
                         int *L12_Jj, int *l3_Jj, int *J12_Jj, int *two_j3_Jj, int *S12_Jj, int *T12_Jj,
                         MKL_INT64 Lmax, MKL_INT64 max_L12, MKL_INT64 max_l3, MKL_INT64 two_J, MKL_INT64 two_T,
                         double *SixJ_array, int two_jmax_SixJ, double* Gtilde_store);

/* Function defined for WPCD method, adapted simply from generate_Ptilde_new */
double calculate_P123_element_in_WP_basis ( int  alpha_idx, int  p_idx_WP, int  q_idx_WP, 
                                            int alphap_idx, int pp_idx_WP, int qp_idx_WP, 
                                            int Np_per_WP, double *p_array, double *wp_array,
                                            int Nq_per_WP, double *q_array, double *wq_array,
                                            int Nx,        double *x_array, double *wx_array,
                                            int Np_WP,     double *p_array_WP_bounds,
                                            int Nq_WP,     double *q_array_WP_bounds,
                                            int Nalpha,
                                            double* Gtilde_store );


double Atilde (int alpha, int alphaprime, int Ltotal, int Jj_dim, int *L12_Jj, int *l3_Jj,
               int *J12_Jj, int *two_j3_Jj, int *S12_Jj, int *T12_Jj, int two_J, int two_T,
               double *SixJ_array, int two_jmax_SixJ);

double Gtilde_new (double p, double q, double x, int alpha, int alphaprime, int N_alpha,
                   int Lmax, int *L12_Jj, int *l3_Jj, double *A_store, int two_Jtotal);

double Gtilde_subarray_new (double p,
							double q,
							double x,
							int L12, int L12prime, int max_L12,
							int l3,  int l3prime,  int max_l3,
							double *A_store_alpha_alphaprime_subarray,
							int two_Jtotal,
							double* ClebschGordan_data,
							double* gsl_Plm_1_subarray,
							double* gsl_Plm_2_subarray,
							double* gsl_Plm_3_subarray,
                            double* prefac_L_array,
                            double* prefac_l_array,
							int two_jmax_Clebsch);



#endif // AUXILIARY_H