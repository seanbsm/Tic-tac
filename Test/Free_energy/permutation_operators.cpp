
#include "permutation_operators.h"

#define FILE "../../../Data/Faddeev_permutation_operator_benchmark/P123_J3_1_PAR_1_T3_1.h5"

void read_P123_data_file(double* P123){

    hid_t       file_id, dataset_id;  /* identifiers */
    herr_t      status;
    int         i, j, dset_data[4][6];
 
    /* Initialize the dataset. */
    for (i = 0; i < 4; i++)
       for (j = 0; j < 6; j++)
          dset_data[i][j] = i * 6 + j + 1;
 
    /* Open an existing file. */
    file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
 
    /* Open an existing dataset. */
    dataset_id = H5Dopen2(file_id, "/dset", H5P_DEFAULT);
 
    /* Write the dataset. */
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                      dset_data);
 
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                     dset_data);
 
    /* Close the dataset. */
    status = H5Dclose(dataset_id);
 
    /* Close the file. */
    status = H5Fclose(file_id);
}

void calculate_antisymmetrization_operator(double* state_3N_bra_array,
                                           double* state_3N_ket_array,
                                           int &Np, double* p_array, double* wp_array,
                                           int &Nq, double* q_array, double* wq_array,
                                           int& Nalpha, int* L_2N, int* S_2N, int* J_2N, int* T_2N, int* l_3N, int* two_j_3N,
                                           double* A123){
    
    double* P123 = NULL;

    read_P123_data_file(P123);


    /*int D123_dim = Np_3N * Nq_3N * Jj_dim;
    int D123_dim_sq = D123_dim * D123_dim;

    float *P123_float;
    P123_float = new float[D123_dim_sq];

    printf("chop matrix elements...\n");


    double eps = 1e-10;
    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
        for (int q_index = 0; q_index <= Nq_3N - 1; q_index++){
            for (int p_index = 0; p_index <= Np_3N - 1; p_index++){
                for (int alpha_prime = 0; alpha_prime <= Jj_dim - 1; alpha_prime++){
                    for (int q_prime_index = 0; q_prime_index <= Nq_3N - 1; q_prime_index++){
                        for (int p_prime_index = 0; p_prime_index <= Np_3N - 1; p_prime_index++){
                            
                            int index = (int) (alpha * Nq_3N * Np_3N + q_index * Np_3N + p_index) * D123_dim +
                                              alpha_prime * Nq_3N * Np_3N + q_prime_index * Np_3N + p_prime_index;

                            if (fabs(P123_store[index]) < eps) {
                                P123_float[index] = 0.0;
                            }
                            else {
                                P123_float[index] = P123_store[index];
                            }
                        }
                    }
                }
            }
        }
    }*/
    
}