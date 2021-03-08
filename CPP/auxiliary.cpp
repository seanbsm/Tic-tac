
#include "auxiliary.h"

double P (double x, void *params)       // P, dP, fdf necessary for initalizing of the object FDF
{
    int n = *(int *) params;
    return (gsl_sf_legendre_Pl(n, x));
}

double dP (double x, void *params)
{
    int n = *(int *) params;
    return ((n * x * gsl_sf_legendre_Pl(n, x) - n * gsl_sf_legendre_Pl(n - 1, x)) / (x * x - 1));
}

void fdf (double x, void *params, double *y, double *dy)
{

    int n = *(int *) params;

    *y = gsl_sf_legendre_Pl(n, x);
    *dy = (n * x * gsl_sf_legendre_Pl(n, x) - n * gsl_sf_legendre_Pl(n - 1, x)) / (x * x - 1);
}

double Legendre (double x, int n)               // P(n,x)
{
    return (gsl_sf_legendre_Pl(n, x));
}

double dLegendre (double x, int n)              // P'(n,x)
{
    return ((n * x * gsl_sf_legendre_Pl(n, x) - n * gsl_sf_legendre_Pl(n - 1, x)) / (x * x - 1));
}

void roots (double *x, int n)                   // calculate roots of legendre polynomials
{

    int status;
    double q[n];

    for (int i = 0; i <= n - 1; i++) q[n - 1 - i] = cos(3.14159 * (i + 1 - 0.25) / (n + 0.5));

    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0;
    gsl_function_fdf FDF;

    for (int i = 0; i <= n - 1; i++)
    {
        FDF.f = &P;
        FDF.df = &dP;
        FDF.fdf = &fdf;
        FDF.params = &n;

        T = gsl_root_fdfsolver_newton;
        s = gsl_root_fdfsolver_alloc (T);
        gsl_root_fdfsolver_set (s, &FDF, q[i]);

        int iter = 0, max_iter = 100;
        do
        {
            iter++;
            status = gsl_root_fdfsolver_iterate (s);
            x0 = q[i];
            q[i] = gsl_root_fdfsolver_root(s);
            status = gsl_root_test_delta (q[i], x0, 0, 1e-9);
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fdfsolver_free(s);
    }
    for (int i = 0; i <= n - 1; i ++)
    {
        x[i] = q[i];
    }
}

void calc_gauss_points (double *q_par, double *w_par, double pmin_par, double pmax_par, int N_par)
{

    double qtemp[N_par];
    roots(qtemp, N_par);

    for (int j = 0; j <= N_par - 1; j++)
    {
        w_par[j] = 2.0 / (1 - qtemp[j] * qtemp[j]) * 1.0 / (dLegendre(qtemp[j], N_par) * dLegendre(qtemp[j], N_par)) * (pmax_par - pmin_par) / 2.0;
        q_par[j] = (pmax_par - pmin_par) * (qtemp[j] + 1) / 2.0 + pmin_par;
    }
}

// transform \int_a^infinity f(x) dx = \int_{arctan(a)}^(Pi/2) f(tan(t)) (1+tan(t)^2) dt
// i.e. transform
void calc_gauss_points_tan (double *q_par, double *w_par, double a, double b, int N_par)
{

    double qtemp[N_par];
    roots(qtemp, N_par);

    double pmin = atan(a);
    double pmax = atan(b);

    for (int j = 0; j <= N_par - 1; j++)
    {
        q_par[j] = tan((pmax - pmin) * (qtemp[j] + 1) / 2.0 + pmin);
        w_par[j] = 2.0 / (1 - qtemp[j] * qtemp[j]) * 1.0 / (dLegendre(qtemp[j], N_par) * dLegendre(qtemp[j], N_par)) * (pmax - pmin) / 2.0 * (1.0 + q_par[j] * q_par[j]);
    }
}

// enforce fewer points around the critical areas around x=1 and x=-1
void calc_gauss_points_ang_opt (double *q_par, double *w_par, double pmin_par, double pmax_par, int N_par)
{
    double qtemp[N_par - 2];
    roots(qtemp, N_par - 2);

    w_par[0] = 2.0 / (1 - 0.0 * 0.0) * 1.0 / (dLegendre(0.0, 1) * dLegendre(0.0, 1)) * (pmin_par + 1) / 2.0;
    q_par[0] = (pmin_par + 1) * (0 + 1) / 2 - 1;

    for (int j = 0; j <= N_par - 3; j++)
    {
        w_par[j + 1] = 2.0 / (1 - qtemp[j] * qtemp[j]) * 1.0 / (dLegendre(qtemp[j], N_par - 2) * dLegendre(qtemp[j], N_par - 2)) * (pmax_par - pmin_par) / 2.0;
        q_par[j + 1] = (pmax_par - pmin_par) * (qtemp[j] + 1) / 2 + pmin_par;
    }

    w_par[N_par - 1] = 2.0 / (1 - 0.0 * 0.0) * 1.0 / (dLegendre(0.0, 1) * dLegendre(0.0, 1)) * (1 - pmax_par) / 2.0;
    q_par[N_par - 1] = (1 - pmax_par) * (0 + 1) / 2 + pmax_par;
}

double ClebschGordan (int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3)
{
    double result = 0.0;
    if ((abs(two_m1) <= two_j1) && (abs(two_m2) <= two_j2) && (abs(two_m1) <= two_j1))
    {
        result = gsl_sf_pow_int(-1, (two_j1 - two_j2 + two_m3) / 2) * sqrt(two_j3 + 1) * gsl_sf_coupling_3j (two_j1, two_j2, two_j3, two_m1, two_m2, -two_m3);
    }
    return result;
}

double ClebschGordan_array (int two_jmax, double *data_array, int two_j1, int two_j2, int two_j3, int two_m1, int two_m2)
{
    // do all index checking while tabulating symbols and restrict sums over M and L so that all conditions are fulfilled
    return data_array[((two_j1 + 1) * (two_j1 + 1) + two_m1 - two_j1 - 1) * (two_jmax + 1) * (two_jmax + 1) * (2 * two_jmax + 1) + ((two_j2 + 1) * (two_j2 + 1) + two_m2 - two_j2 - 1) * (2 * two_jmax + 1) + two_j3];
}

double ThreeJarray (int jmax, double *data_array, int j1, int j2, int j3, int m1, int m2)
{
    // do all index checking while tabulating symbols
    return data_array[((j1 + 1) * (j1 + 1) + m1 - j1 - 1) * (jmax + 1) * (jmax + 1) * (jmax + 1) + ((j2 + 1) * (j2 + 1) + m2 - j2 - 1) * (jmax + 1) + j3];
}

double SixJSymbol (double *SixJ_array, int two_jmax, int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6)
{
    return SixJ_array[
               two_j1 * (two_jmax + 1) * (two_jmax + 1) * (two_jmax + 1) * (two_jmax + 1) * (two_jmax + 1)
               + two_j2 * (two_jmax + 1) * (two_jmax + 1) * (two_jmax + 1) * (two_jmax + 1)
               + two_j3 * (two_jmax + 1) * (two_jmax + 1) * (two_jmax + 1)
               + two_j4 * (two_jmax + 1) * (two_jmax + 1)
               + two_j5 * (two_jmax + 1)
               + two_j6];
}


double NineJSymbol (int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6, int two_j7, int two_j8)
{
    return gsl_sf_coupling_9j(two_j1, two_j2, two_j3, two_j4, 1, two_j5, two_j6, two_j7, two_j8);
}

double ClebschGordan_array (int two_jmax, double *data_array, int two_j1, int two_j2, int two_j3, int two_m1, int two_m2);

int KroneckerDelta(int i, int j)
{
    return (i == j);
}

// actually not needed....
complex<double> Ylm (int l, int m, double theta, double phi)
{

    // inefficient, actually, we only need Ylm at phi=0, associated LegendrePolynomials are sufficient
    if (abs(m) > l) cout << "Ylm error! abs(m) > l\n";

    gsl_complex phi_exp, result;
    complex <double> result_cast;
    int m_abs = abs(m);

    GSL_SET_COMPLEX(&phi_exp, 0.0, m_abs * phi);
    phi_exp = gsl_complex_exp(phi_exp);

    result = gsl_complex_mul_real (phi_exp, gsl_sf_legendre_sphPlm (l, m_abs, cos(theta)));
    if (m < 0) result = gsl_complex_mul_real(gsl_complex_conjugate(result), gsl_sf_pow_int(-1, m));

    result_cast = complex <double> (GSL_REAL(result), GSL_IMAG(result));

    return result_cast;

}

double Plm (int l, int m, double costheta)   // associate Legendre Polynomials with proper normalization, so Ylm = e^{i m phi} Plm
{

    double result = 0;
    if (m >= 0) result = gsl_sf_legendre_sphPlm(l, m, costheta);
    else
    {
        result = gsl_sf_pow_int(-1.0, m) * gsl_sf_legendre_sphPlm(l, -m, costheta);
    }
    return result;
}

// calc interpolating function: p_par(p_par_index) < x < p_par(p_par_index+1)
// f(x) = \sum_i f(p_par(i)) S_fct(i,x,p_par,N_par,S_vec)
double S_fct(int i, double x, double *p_par, int N_par, double *S_vec_par)
{

    int j_index = 0; while ((p_par[j_index] < x) && (j_index < N_par - 1)) j_index++; if (j_index > 0) j_index--;
    double ret;

    if ((x < p_par[N_par - 1] + 1e-5))
    {
        ret = (S_vec_par[j_index * N_par + i]
               + (x - p_par[j_index]) * S_vec_par[j_index * N_par + i + N_par * (N_par - 1)]
               + (x - p_par[j_index]) * (x - p_par[j_index]) * S_vec_par[j_index * N_par + i + 2 * N_par * (N_par - 1)]
               + (x - p_par[j_index]) * (x - p_par[j_index]) * (x - p_par[j_index]) * S_vec_par[j_index * N_par + i + 3 * N_par * (N_par - 1)]
              );
    }
    else ret = 0.0;

    return ret;
}

// more efficient method, assumes that j_index is already known, saves the internal loop
double S_fct(int i, double x, int j_index, double *p_par, int N_par, double *S_vec_par)
{
    double ret;

    if ((x < p_par[N_par - 1] + 1e-5))
    {
        ret = (S_vec_par[j_index * N_par + i]
               + (x - p_par[j_index]) * S_vec_par[j_index * N_par + i + N_par * (N_par - 1)]
               + (x - p_par[j_index]) * (x - p_par[j_index]) * S_vec_par[j_index * N_par + i + 2 * N_par * (N_par - 1)]
               + (x - p_par[j_index]) * (x - p_par[j_index]) * (x - p_par[j_index]) * S_vec_par[j_index * N_par + i + 3 * N_par * (N_par - 1)]
              );
    }
    else ret = 0.0;

    return ret;
}

// cubic spline
void S_spline (double *S, double *p_par, int N_par)
{

    double h[N_par];
    h[0] = 0.0;
    for (int j = 1; j <= N_par - 1; j++)
    {
        h[j] = p_par[j] - p_par[j - 1];
    }

    double lambda[N_par - 1];
    double mu[N_par - 1];
    for (int j = 0; j <= N_par - 2; j++)
    {
        lambda[j] = h[j + 1] / (h[j] + h[j + 1]);
        mu[j] = 1 - lambda[j];
    }

    double p[N_par - 1];
    double q[N_par - 1];

    q[0] = 0.0;
    for (int j = 1; j <= N_par - 2; j++)
    {
        p[j] = mu[j] * q[j - 1] + 2.0;
        q[j] = -lambda[j] / p[j];
    }

    double B[N_par * (N_par - 1)];
    for (int i = 0; i <= N_par - 1; i++)
    {
        for (int j = 0; j <= N_par - 2; j++)
        {
            B[j * N_par + i] =
                KroneckerDelta(i, j + 1) * 6.0 / ((h[j] + h[j + 1]) * h[j + 1])
                - KroneckerDelta(i, j) * 6.0 / (h[j] * h[j + 1])
                + KroneckerDelta(i, j - 1) * 6.0 / ((h[j] + h[j + 1]) * h[j]);
        }
    }

    double A[N_par * (N_par - 1)];
    for (int i = 0; i <= N_par - 1; i++)
    {
        A[0 * N_par + i] = 0.0;
        for (int j = 1; j <= N_par - 2; j++)
        {
            A[j * N_par + i] = B[j * N_par + i] / p[j] - mu[j] / p[j] * A[(j - 1) * N_par + i];
        }
    }

    double C[N_par * N_par];
    for (int i = 0; i <= N_par - 1; i++)
    {
        C[(N_par - 1)*N_par + i] = 0.0;
        for (int j = N_par - 2; j >= 0; j--)
        {
            C[j * N_par + i] = q[j] * C[(j + 1) * N_par + i] + A[j * N_par + i];
        }
    }

    for (int i = 0; i <= N_par - 1; i++)
    {
        for (int j = 0; j <= N_par - 2; j++)
        {
            S[j * N_par + i                          ] = KroneckerDelta(i, j);
            S[j * N_par + i +     N_par * (N_par - 1)] = (KroneckerDelta(i, j + 1) - KroneckerDelta(i, j)) / h[j + 1] - h[j + 1] / 6.0 * (2.0 * C[j * N_par + i] + C[(j + 1) * N_par + i]);
            S[j * N_par + i + 2 * N_par * (N_par - 1)] = 0.5 * C[j * N_par + i];
            S[j * N_par + i + 3 * N_par * (N_par - 1)] = (C[(j + 1) * N_par + i] - C[j * N_par + i]) / (6.0 * h[j + 1]);
        }
    }
}

double lin_poly_interpolation (double x1, double y1, double x2, double y2, double x)
{   
    return ((x - x2) / (x1 - x2) * y1 + (x - x1) / (x2 - x1) * y2);
}

complex<double> clin_poly_interpolation (double x1, complex<double> y1, double x2, complex<double> y2, double x)
{   
    return ((x - x2) / (x1 - x2) * y1 + (x - x1) / (x2 - x1) * y2);
}

double quadr_poly_interpolation (double x1, double y1, double x2, double y2, double x3, double y3, double x)
{   
    return ((x - x2) * (x - x3) / (x1 - x2) / (x1 - x3) * y1 + (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3) * y2 + (x - x1) * (x - x2) / (x3 - x1) / (x3 - x2) * y3);
}

// f1, f2, f3, f4 function values on a rectangle counted counterclockwise beginning in the lower left corner
double lin_poly_2dim (double x1, double x2, double y1, double y2, double f1, double f2, double f3, double f4, double x, double y, double cutoff)
{
    
    double t = (x - x1) / (x2 - x1);
    double u = (y - y1) / (y2 - y1);
    
    if ((x < cutoff) && (y < cutoff))
    {   
        return  (1.0 - t) * (1.0 - u) * f1 + t * (1.0 - u) * f2 + t * u * f3 + (1.0 - t) * u * f4;
    }
    else return 0.0;
}

// f1, f2, f3, f4, f5, f6, f7, f8 function values on a rectangle counted counterclockwise beginning in the lower left corner, first rectangle f1-f4 lower coordinate in z, f5-f8 higher coordinate
double lin_poly_3dim (double x1, double x2, double y1, double y2, double z1, double z2,
                      double f1, double f2, double f3, double f4, double f5, double f6, double f7, double f8,
                      double x, double y, double z, double cutoff)
{

    double xd = (x - x1) / (x2 - x1);
    double yd = (y - y1) / (y2 - y1);
    double zd = (z - z1) / (z2 - z1);

    double i1 = f1 * (1 - zd) + f5 * zd;
    double i2 = f4 * (1 - zd) + f8 * zd;
    double j1 = f2 * (1 - zd) + f6 * zd;
    double j2 = f3 * (1 - zd) + f7 * zd;

    double w1 = i1 * (1 - yd) + i2 * yd;
    double w2 = j1 * (1 - yd) + j2 * yd;

    if ((x < cutoff) && (y < cutoff))
    {
        return  w1 * (1 - xd) + w2 * xd;
    }
    else return 0.0;
}

double Floc_l (int l1, int l2, int l1prime, int l2prime, int l, double p, double q, double pprime, double qprime, int Npq, double *ClebschGordan_data, int two_jmax, double *p_tilde_vec, double *wp_tilde_vec, double *q_tilde_vec, double *wq_tilde_vec, double *F_interpolate, int lmax, double *Plm_p_1, double *Plm_p_2, double *Plm_q_1, double *Plm_q_2)
{

    double result = 0.0;
    double prefac = 1.0 / (p * pprime * q * qprime);

    int two_l1 = 2 * l1;
    int two_l2 = 2 * l2;
    int two_l1prime = 2 * l1prime;
    int two_l2prime = 2 * l2prime;
    int two_l = 2 * l;

    // tabulate vector spherical harmonics at phi=0 (real quantity)
    double Yllm_p[Npq];
    double Yllm_q[Npq];

    for (int i = 0; i <= Npq - 1; i++)
    {

        Yllm_p[i] = 0.0;
        for (int ml = -GSL_MIN_INT(l1, l1prime); ml <= GSL_MIN_INT(l1, l1prime); ml++)
        {
            Yllm_p[i] += ClebschGordan_array(two_jmax, ClebschGordan_data, two_l1prime, two_l1, two_l, 2 * ml, -2 * ml) * Plm_p_2[((l1prime + 1) * (l1prime + 1) + ml - l1prime - 1) * Npq + i] * Plm_p_1[((l1 + 1) * (l1 + 1) - ml - l1 - 1) * Npq + i];
        }
    }

    for (int j = 0; j <= Npq - 1; j++)
    {

        Yllm_q[j] = 0.0;
        for (int ml = -GSL_MIN_INT(l2, l2prime); ml <= GSL_MIN_INT(l2, l2prime); ml++)
        {
            Yllm_q[j] += ClebschGordan_array(two_jmax, ClebschGordan_data, two_l2prime, two_l2, two_l, 2 * ml, -2 * ml) * Plm_q_2[((l2prime + 1) * (l2prime + 1) + ml - l2prime - 1) * Npq + j] * Plm_q_1[((l2 + 1) * (l2 + 1) - ml - l2 - 1) * Npq + j];
        }
    }

    for (int i = 0; i <= Npq - 1; i++)
    {
        for (int j = 0; j <= Npq - 1; j++)
        {
            result += wp_tilde_vec[i] * wq_tilde_vec[j] * p_tilde_vec[i] * q_tilde_vec[j]
                      * Yllm_p[i] * Yllm_q[j]
                      * F_interpolate[l * Npq * Npq + i * Npq + j]
                      ;
        } // j
    } // i

    return prefac * result;

}

double Gtilde (double p, double q, double x, int alpha, int alphaprime, int Jj_dim, int Lmax, int kLegendremax, double *g_array_par, int *L12_Jj, int *l3_Jj)
{

    double sum = 0.0;

    //cout << "kLegenedremax = " << kLegendremax << "\n";

    for (int k = 0; k <= kLegendremax; k++)
    {
        for (int L1 = 0; L1 <= L12_Jj[alphaprime]; L1++)
        {
            for (int l1 = 0; l1 <= l3_Jj[alphaprime]; l1++)
            {
                sum += gsl_sf_legendre_Pl(k, x)
                       * pow(p, L1 + l1)
                       * pow(q, L12_Jj[alphaprime] - L1 + l3_Jj[alphaprime] - l1)
                       * g_array_par[
                           k * (Lmax + 1) * (Lmax + 1) * Jj_dim * Jj_dim
                           + L1 * (Lmax + 1) * Jj_dim * Jj_dim
                           + l1 * Jj_dim * Jj_dim
                           + alpha * Jj_dim
                           + alphaprime];
            }
        }
    }
    return sum;
}

double pi1_tilde(double p, double q, double x)
{
    return sqrt(0.25 * p * p + 9.0 / 16.0 * q * q + 0.75 * p * q * x);
}


double pi2_tilde(double p, double q, double x)
{
    return sqrt(p * p + 0.25 * q * q - p * q * x);
}

double pi1_prime_tilde(double p, double q, double x)
{
    return sqrt(0.25 * p * p + 9.0 / 16.0 * q * q - 0.75 * p * q * x);
}


double pi2_prime_tilde(double p, double q, double x)
{
    return sqrt(p * p + 0.25 * q * q + p * q * x);
}

void generate_Ptilde_new (double *P123_store, MKL_INT64 Pdim, MKL_INT64 N_p, double *p, MKL_INT64 N_q, double *q, MKL_INT64 N_x, double *x, double *wx, MKL_INT64 Jj_dim, double pmax, double qmax, int *L12_Jj, int *l3_Jj, int *J12_Jj, int *two_j3_Jj, int *S12_Jj, int *T12_Jj, MKL_INT64 Lmax, MKL_INT64 max_L12, MKL_INT64 max_l3, MKL_INT64 two_J, MKL_INT64 two_T, double *SixJ_array, int two_jmax_SixJ/*, double *Atilde_store, double *Btilde_store*/, double* Gtilde_store)
{

    // unused variable: int iZERO = 0;

    for (MKL_INT64 index = 0; index <= Pdim * Pdim - 1; index++)
    {
        P123_store[index] = 0.0;
    }

    double S_vec_p[4 * N_p * (N_p - 1)];
    S_spline(S_vec_p, p, N_p);

    double S_vec_q[4 * N_q * (N_q - 1)];
    S_spline(S_vec_q, q, N_q);

    double interpolate_p;
    double interpolate_q;

    // unused variable: int bary_d = 1;

    // unused variable: double bary_wp[N_p];
    // unused variable: double bary_wq[N_q];

    #pragma omp parallel private(interpolate_p, interpolate_q)
    {
        #pragma omp for

        for (MKL_INT64 c = 0; c <= Pdim - 1; c++)
        {

            int c_fortran = c + 1;
            int c_global = c_fortran - 1; //indxl2g_(&c_fortran, &X123_Nb, &mycol123, &iZERO, &npcols123)-1;
            // unused variable: int c_store = c_global;

            int alphaprime = (int) (c_global / (N_q * N_p));
            c_global = c_global - alphaprime * N_q * N_p;

            int qprime_index = (int) (c_global / (N_p));
            c_global = c_global - qprime_index * N_p;

            int pprime_index = c_global;

	    /*
            if ((alphaprime >= Jj_dim) || (alphaprime < 0))
            {
                cout << "alphaprime problem!\n";
                cout << alphaprime << "\n";
                //MPI_Finalize();
            }

            if ((qprime_index >= N_q) || (qprime_index < 0))
            {
                cout << "qprime problem\n";
                cout << qprime_index << "\n";
                //MPI_Finalize();
            }

            if ((pprime_index >= N_p) || (pprime_index < 0))
            {
                cout << "pprime problem\n";
                //MPI_Finalize();
            }
	    */

            for (MKL_INT64 r = 0; r <= Pdim - 1; r++)
            {

                int r_fortran = r + 1;
                int r_global = r_fortran - 1;//indxl2g_(&r_fortran, &X123_Mb, &myrow123, &iZERO, &nprows123)-1;
                // unused variable: int r_store = r_global;

                int alpha = (int) (r_global / (N_q * N_p));
                r_global = r_global - alpha * N_q * N_p;

                int q_index = (int) (r_global / (N_p));
                r_global = r_global - q_index * N_p;

                int p_index = r_global;

		/*
                if ((alpha >= Jj_dim) || (alpha < 0))
                {
                    cout << "alpha problem!\n";
                    //MPI_Finalize();
                }

                if ((q_index >= N_q) || (q_index < 0))
                {
                    cout << "q problem\n";
                    //MPI_Finalize();
                }

                if ((p_index >= N_p) || (p_index < 0))
                {
                    cout << "p problem\n";
                    //MPI_Finalize();
                }
		*/

                for (int x_index = 0; x_index <= N_x - 1; x_index++)
                {

                    double pi1 = pi1_tilde(p[p_index], q[q_index], x[x_index]);
                    double pi2 = pi2_tilde(p[p_index], q[q_index], x[x_index]);

                    int pi1_index = 0; while ((p[pi1_index] < pi1) && (pi1_index < N_p - 1)) pi1_index++; if (pi1_index > 0) pi1_index--;
                    int pi2_index = 0; while ((q[pi2_index] < pi2) && (pi2_index < N_q - 1)) pi2_index++; if (pi2_index > 0) pi2_index--;

                    double costheta1 = -(0.5 * p[p_index] + 0.75 * q[q_index] * x[x_index]) / pi1;
                    double costheta2 = (p[p_index] - 0.5 * q[q_index] * x[x_index]) / pi2;

                    if (fabs(costheta1) > 1) cout << "costheta1 problem: " << costheta1 << "\n";
                    if (fabs(costheta2) > 1) cout << "costheta2 problem: " << costheta2 << "\n";

                    if (pi1 <= p[N_p - 1]) // avoid artefacts at the upper limits of mesh system
                    {

                        // spline

                        interpolate_p = (
                                            S_vec_p[pi1_index * N_p + pprime_index]
                                            + (pi1 - p[pi1_index]) * S_vec_p[pi1_index * N_p + pprime_index + N_p * (N_p - 1)]
                                            + (pi1 - p[pi1_index]) * (pi1 - p[pi1_index]) * S_vec_p[pi1_index * N_p + pprime_index + 2 * N_p * (N_p - 1)]
                                            + (pi1 - p[pi1_index]) * (pi1 - p[pi1_index]) * (pi1 - p[pi1_index]) * S_vec_p[pi1_index * N_p + pprime_index + 3 * N_p * (N_p - 1)]
                                        );

                    }
                    else
                    {
                        interpolate_p = 0.0;
                    }

                    if (pi2 <= q[N_q - 1]) // avoid spline artefacts at the upper limits of mesh system
                    {

                        // spline

                        interpolate_q = (
                                            S_vec_q[pi2_index * N_q + qprime_index]
                                            + (pi2 - q[pi2_index]) * S_vec_q[pi2_index * N_q + qprime_index + N_q * (N_q - 1)]
                                            + (pi2 - q[pi2_index]) * (pi2 - q[pi2_index]) * S_vec_q[pi2_index * N_q + qprime_index + 2 * N_q * (N_q - 1)]
                                            + (pi2 - q[pi2_index]) * (pi2 - q[pi2_index]) * (pi2 - q[pi2_index]) * S_vec_q[pi2_index * N_q + qprime_index + 3 * N_q * (N_q - 1)]
                                        );

                    }
                    else
                    {
                        interpolate_q = 0.0;
                    }

                    // new:
                    P123_store[r * Pdim + c] += wx[x_index] * interpolate_p * interpolate_q * Gtilde_store[alpha * Jj_dim * N_p * N_q * N_x + alphaprime * N_p * N_q * N_x + p_index * N_q * N_x + q_index * N_x + x_index];

                    // old:
                    /*
                    // note fortran convention for rows/columns! needed to make descriptors work properly
                    for (MKL_INT64 Ltotal = max(abs((int) (L12_Jj[alpha] - l3_Jj[alpha])), abs((int) (L12_Jj[alphaprime] - l3_Jj[alphaprime]))); Ltotal <= min((int) ((two_J + 5) / 2), min((int) (L12_Jj[alpha] + l3_Jj[alpha]), (int) (L12_Jj[alphaprime] + l3_Jj[alphaprime]))); Ltotal++)
                    {

                        MKL_INT64 index = L12_Jj[alpha] * (max_l3 + 1) * (max_L12 + 1) * (max_l3 + 1) * N_p * N_q * N_x * (Lmax + 1)
                                          + l3_Jj[alpha] * (max_L12 + 1) * (max_l3 + 1) * N_p * N_q * N_x * (Lmax + 1)
                                          + L12_Jj[alphaprime] * (max_l3 + 1) * N_p * N_q * N_x * (Lmax + 1)
                                          + l3_Jj[alphaprime] * N_p * N_q * N_x * (Lmax + 1)
                                          + p_index * N_q * N_x * (Lmax + 1)
                                          + q_index * N_x * (Lmax + 1)
                                          + x_index * (Lmax + 1)
                                          + Ltotal;

                        if (index < 0) cout << "WTF! : " << index << " " << L12_Jj[alpha] << " " << (max_l3 + 1) << " " << (max_L12 + 1) << " " << N_p << " " << N_q << " " << N_x << "\n";

                        P123_store[r * Pdim + c] += // note fortran convention for rows/columns! needed to make descriptors work properly

                            Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal]
                            * 8.0 * M_PI * M_PI * wx[x_index] * interpolate_p * interpolate_q
                            * Btilde_store[index]
                            ;
                    } // Ltotal
                    */
                } // x_index
            } // j
        } // i
    } // pragma
}

double calculate_P123_element_in_WP_basis ( int  alpha_idx, int  p_idx_WP, int  q_idx_WP, 
                                            int alphap_idx, int pp_idx_WP, int qp_idx_WP, 
                                            int Np_per_WP, double *p_array, double *wp_array,
                                            int Nq_per_WP, double *q_array, double *wq_array,
                                            int Nx,        double *x_array, double *wx_array,
                                            int Np_WP,     double *p_array_WP_bounds,
                                            int Nq_WP,     double *q_array_WP_bounds,
                                            int Nalpha,
                                            double* Gtilde_store ){

    bool   run_tests    = false;
    double integral_sum = 0;

    double WP_q_bound_lower = q_array_WP_bounds[q_idx_WP];
    double WP_q_bound_upper = q_array_WP_bounds[q_idx_WP+1];
    double WP_p_bound_lower = p_array_WP_bounds[p_idx_WP];
    double WP_p_bound_upper = p_array_WP_bounds[p_idx_WP+1];

    int pp_idx_lower = Np_per_WP* pp_idx_WP;
    int pp_idx_upper = Np_per_WP*(pp_idx_WP+1);

    int qp_idx_lower = Nq_per_WP* qp_idx_WP;
    int qp_idx_upper = Nq_per_WP*(qp_idx_WP+1);

    int Np = Np_per_WP*Np_WP;
    int Nq = Nq_per_WP*Nq_WP;
    
    double pp = 0;
    double qp = 0;
    double wp_p_f = 0;
    double wq_q_f = 0;
    double integral_factors = 0;
    double p_bar = 0;
    double q_bar = 0;
    double x     = 0;
    double G_element = 0;

    int G_pp_idx = 0;
    int G_pp_qp_idx = 0;

    for (int pp_idx=pp_idx_lower; pp_idx<pp_idx_upper; pp_idx++){
        /* p' momentum */
        pp = p_array[pp_idx];
        /* Numerical integral coefficients */
        wp_p_f = wp_array[pp_idx]*pp*p_weight_function(pp);

        G_pp_idx = (pp_idx-pp_idx_lower)*Nq_per_WP*Nx;

        for (int qp_idx=qp_idx_lower; qp_idx<qp_idx_upper; qp_idx++){
            /* q' momentum */
            qp = q_array[qp_idx];
            /* Numerical integral coefficients */
            wq_q_f = wq_array[qp_idx]*qp*q_weight_function(qp);

            G_pp_qp_idx = G_pp_idx + (qp_idx-qp_idx_lower)*Nx;

            integral_factors = wp_p_f * wq_q_f;

            for (int x_idx=0; x_idx<Nx; x_idx++){
                x = x_array[x_idx];

                p_bar = pi1_tilde(pp, qp, x);
                /* Check if p_bar is in bin p_idx_WP, if not then move on to next loop iteration */
                if ( p_bar<WP_p_bound_lower or WP_p_bound_upper<p_bar ){
                    continue;
                }

                q_bar = pi2_tilde(pp, qp, x);
                /* Check if q_bar is in bin q_idx_WP, if not then move on to next loop iteration */
                if ( q_bar<WP_q_bound_lower or WP_q_bound_upper<q_bar ){
                    continue;
                }
                
                if (run_tests){
                    double costheta1 = -(0.5 * pp + 0.75 * qp * x) / p_bar;
                    double costheta2 =        (pp - 0.5  * qp * x) / q_bar;
                    if (fabs(costheta1) > 1) cout << "costheta1 problem: " << costheta1 << "\n";
                    if (fabs(costheta2) > 1) cout << "costheta2 problem: " << costheta2 << "\n";
                }

                //G_element = Gtilde_store[alphap_idx*Nalpha*Np*Nq*Nx + alpha_idx*Np*Nq*Nx + pp_idx*Nq*Nx + qp_idx*Nx + x_idx];
                G_element = Gtilde_store[G_pp_qp_idx + x_idx];
                
                /* BEWARE: I'M NOT 100% SURE ON THE GTILDE ALPHA-INDEXING */
                integral_sum += integral_factors * wx_array[x_idx] * G_element / (p_bar*q_bar);
            } // x_idx
        } // qp_idx
    } // pp_idx

    double WP_qp_bound_lower = q_array_WP_bounds[qp_idx_WP];
    double WP_qp_bound_upper = q_array_WP_bounds[qp_idx_WP+1];
    double WP_pp_bound_lower = p_array_WP_bounds[pp_idx_WP];
    double WP_pp_bound_upper = p_array_WP_bounds[pp_idx_WP+1];

    double norm_WP = q_normalization( WP_q_bound_lower,  WP_q_bound_upper)
                    *p_normalization( WP_p_bound_lower,  WP_p_bound_upper)
                    *q_normalization(WP_qp_bound_lower, WP_qp_bound_upper)
                    *p_normalization(WP_pp_bound_lower, WP_pp_bound_upper);

    return integral_sum / norm_WP;
}

double Atilde (int alpha, int alphaprime, int Ltotal, int Jj_dim, int *L12_Jj, int *l3_Jj, int *J12_Jj, int *two_j3_Jj, int *S12_Jj, int *T12_Jj, int two_J, int two_T, double *SixJ_array, int two_jmax_SixJ)
{

    double ret = 0.0;

    for (int two_Stotal = 1; two_Stotal <= 3; two_Stotal += 2)
    {

        ret +=          sqrt(    (2 * J12_Jj[alpha] + 1) *     (two_j3_Jj[alpha] + 1) *     (2 * S12_Jj[alpha] + 1) *     (2 * T12_Jj[alpha] + 1)
                                 * (2 * J12_Jj[alphaprime] + 1) * (two_j3_Jj[alphaprime] + 1) * (2 * S12_Jj[alphaprime] + 1) * (2 * T12_Jj[alphaprime] + 1)
                            )
                        * gsl_sf_pow_int(-1, S12_Jj[alphaprime] + T12_Jj[alphaprime])
                        * (two_Stotal + 1)
                        * gsl_sf_coupling_9j(2 * L12_Jj[alpha],     2 * S12_Jj[alpha],     2 * J12_Jj[alpha],     2 * l3_Jj[alpha],     1, two_j3_Jj[alpha],     2 * Ltotal, two_Stotal, two_J)
                        * gsl_sf_coupling_9j(2 * L12_Jj[alphaprime], 2 * S12_Jj[alphaprime], 2 * J12_Jj[alphaprime], 2 * l3_Jj[alphaprime], 1, two_j3_Jj[alphaprime], 2 * Ltotal, two_Stotal, two_J)
                        * SixJSymbol(SixJ_array, two_jmax_SixJ, 1, 1, 2 * S12_Jj[alpha], 1, two_Stotal, 2 * S12_Jj[alphaprime])
                        * SixJSymbol(SixJ_array, two_jmax_SixJ, 1, 1, 2 * T12_Jj[alpha], 1, two_T, 2 * T12_Jj[alphaprime]);

    }

    return ret;
}

double Gtilde_new (double p, double q, double x, int alpha, int alphaprime, int N_alpha, int Lmax, int *L12_Jj, int *l3_Jj, double *A_store, int two_Jtotal)
{

    double ret = 0.0;

    double fac1, fac2;

    double pi1 = pi1_tilde(p, q, x);
    double pi2 = pi2_tilde(p, q, x);

    double costheta1 = -(0.5 * p + 0.75 * q * x) / pi1;
    double costheta2 = (p - 0.5 * q * x) / pi2;

    /* Prevent numerical error in Plm */
    if ( costheta1>1 ){
        costheta1 = 1;
    }
    else if( costheta1<-1 ){
        costheta1 = -1;
    }
    
    /* Prevent numerical error in Plm */
    if ( costheta2>1 ){
        costheta2 = 1;
    }
    else if( costheta2<-1 ){
        costheta2 = -1;
    }

    int L12 = L12_Jj[alpha];
    int l3 = l3_Jj[alpha];
    int L12prime = L12_Jj[alphaprime];
    int l3prime = l3_Jj[alphaprime];

    for (int Ltotal = max(abs(L12 - l3), abs(L12prime - l3prime)); Ltotal <= min((two_Jtotal + 5) / 2, min(L12 + l3, L12prime + l3prime)); Ltotal++)
    {
        fac1 = 8.0 * M_PI * M_PI * A_store[alpha * N_alpha * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal];
        
        for (int Mtotal = -min(l3, Ltotal); Mtotal <= min(l3, Ltotal); Mtotal++)
        {

            fac2 = ClebschGordan(2 * L12, 2 * l3, 2 * Ltotal, 0, 2 * Mtotal, 2 * Mtotal)
                   * sqrt((2.0 * L12 + 1) / (4 * M_PI))
                   * gsl_sf_pow_int(-1, Mtotal)
                   * Plm(l3, Mtotal, x); // -1^M phase since azimutal angles of p' and q' = pi

            for (int M12primesum = -L12prime; M12primesum <= L12prime; M12primesum++)
            {
                if (abs(Mtotal - M12primesum) <= l3prime)
                {
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

double Gtilde_subarray_new (double p, double q, double x, int L12, int L12prime, int l3, int l3prime, double *A_store_alpha_alphaprime_subarray, int two_Jtotal){

    double ret = 0.0;

    double fac1, fac2;

    double pi1 = pi1_tilde(p, q, x);
    double pi2 = pi2_tilde(p, q, x);

    double costheta1 = -(0.5 * p + 0.75 * q * x) / pi1;
    double costheta2 = (p - 0.5 * q * x) / pi2;

    /* Prevent numerical error in Plm */
    if ( costheta1>1 ){
        costheta1 = 1;
    }
    else if( costheta1<-1 ){
        costheta1 = -1;
    }
    
    /* Prevent numerical error in Plm */
    if ( costheta2>1 ){
        costheta2 = 1;
    }
    else if( costheta2<-1 ){
        costheta2 = -1;
    }

    for (int Ltotal = max(abs(L12 - l3), abs(L12prime - l3prime)); Ltotal <= min((two_Jtotal + 5) / 2, min(L12 + l3, L12prime + l3prime)); Ltotal++){
        fac1 = 8.0 * M_PI * M_PI * A_store_alpha_alphaprime_subarray[Ltotal];

        for (int Mtotal = -min(l3, Ltotal); Mtotal <= min(l3, Ltotal); Mtotal++)
        {

            fac2 = ClebschGordan(2 * L12, 2 * l3, 2 * Ltotal, 0, 2 * Mtotal, 2 * Mtotal)
                   * sqrt((2.0 * L12 + 1) / (4 * M_PI))
                   * gsl_sf_pow_int(-1, Mtotal)
                   * Plm(l3, Mtotal, x); // -1^M phase since azimutal angles of p' and q' = pi

            for (int M12primesum = -L12prime; M12primesum <= L12prime; M12primesum++)
            {
                if (abs(Mtotal - M12primesum) <= l3prime)
                {
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
