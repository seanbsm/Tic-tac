
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

int kronecker_delta(int i, int j){
    return (i == j);
}

/* Cubic spline */
void S_spline (double* S, double* p_par, int N_par){
    double h[N_par];
    h[0] = 0.0;
    for (int j = 1; j <= N_par - 1; j++){
        h[j] = p_par[j] - p_par[j - 1];
    }

    double lambda[N_par - 1];
    double mu[N_par - 1];
    for (int j = 0; j <= N_par - 2; j++){
        lambda[j] = h[j + 1] / (h[j] + h[j + 1]);
        mu[j] = 1 - lambda[j];
    }

    double p[N_par - 1];
    double q[N_par - 1];

    q[0] = 0.0;
    for (int j = 1; j <= N_par - 2; j++){
        p[j] = mu[j] * q[j - 1] + 2.0;
        q[j] = -lambda[j] / p[j];
    }

    double B[N_par * (N_par - 1)];
    for (int i = 0; i <= N_par - 1; i++){
        for (int j = 0; j <= N_par - 2; j++){
            B[j * N_par + i] =
                  kronecker_delta(i, j + 1) * 6.0 / ((h[j] + h[j + 1]) * h[j + 1])
                - kronecker_delta(i, j) * 6.0 / (h[j] * h[j + 1])
                + kronecker_delta(i, j - 1) * 6.0 / ((h[j] + h[j + 1]) * h[j]);
        }
    }

    double A[N_par * (N_par - 1)];
    for (int i = 0; i <= N_par - 1; i++){
        A[0 * N_par + i] = 0.0;
        for (int j = 1; j <= N_par - 2; j++){
            A[j * N_par + i] = B[j * N_par + i] / p[j] - mu[j] / p[j] * A[(j - 1) * N_par + i];
        }
    }

    double C[N_par * N_par];
    for (int i = 0; i <= N_par - 1; i++){
        C[(N_par - 1)*N_par + i] = 0.0;
        for (int j = N_par - 2; j >= 0; j--){
            C[j * N_par + i] = q[j] * C[(j + 1) * N_par + i] + A[j * N_par + i];
        }
    }

    for (int i = 0; i <= N_par - 1; i++){
        for (int j = 0; j <= N_par - 2; j++){
            S[j * N_par + i                  ] = kronecker_delta(i, j);
            S[j * N_par + i +  N_par * (N_par - 1)] = (kronecker_delta(i, j + 1) - kronecker_delta(i, j)) / h[j + 1] - h[j + 1] / 6.0 * (2.0 * C[j * N_par + i] + C[(j + 1) * N_par + i]);
            S[j * N_par + i + 2 * N_par * (N_par - 1)] = 0.5 * C[j * N_par + i];
            S[j * N_par + i + 3 * N_par * (N_par - 1)] = (C[(j + 1) * N_par + i] - C[j * N_par + i]) / (6.0 * h[j + 1]);
        }
    }
}

void iterate_faddeev(double* state_3N_symm_array,
                     double* P123_array,
                     int Np, double* p_array, double* wp_array,
                     int Nq, double* q_array, double* wq_array,
                     int Nalpha, int* L_2N_array, int* S_2N_array, int* J_2N_array, int* T_2N_array, int* l_3N_array, int* two_j_3N_array,
                     int two_T, int two_J, int PAR){
    
    /* Make q and p integral meshes */

    /* Construct spline arrays */
    double S_vec_p[4*Np*(Np-1)];
    S_spline(S_vec_p, p_array, Np);

    double S_vec_q[4*Nq*(Nq-1)];
    S_spline(S_vec_q, q_array, Nq);

    double interpolate_p;
    double interpolate_q;

    int bary_d = 1;

    double bary_wp[Np];
    double bary_wq[Nq];

    /* First faddeev state summation start */

    /* Loop over alpha prime */
    for (int idx_alpha_p=0; idx_alpha_p<Nalpha; idx_alpha_p++){
        /* Loop over alpha double-prime */
        for (int idx_alpha_pp=0; idx_alpha_pp<Nalpha; idx_alpha_pp++){
            /* Loop over q prime */
            for (int idx_q_p=0; idx_q_p<Nq; idx_q_p++){
                /* Loop over x */
                for (int idx_x=0; idx_x<Nx; idx_x++){
                    
                }
            }
        }
    }

    /* First faddeev state summation end */

}