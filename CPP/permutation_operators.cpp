
#include "permutation_operators.h"

void check_h5_read_call(herr_t ret){
    if (ret<0){
        raise_error("Erroneous control return-value from calling H5Dread()");
    }
}

void check_h5_read_table_call(herr_t ret){
    if (ret<0){
        raise_error("Erroneous control return-value from calling H5TBread_table()");
    }
}

void check_h5_close_call(herr_t ret){
    if (ret<0){
        raise_error("Erroneous control return-value from calling H5Dclose()");
    }
}

void read_P123_h5_data_file(std::string file_path, double* P123, int Nq_3N, double *q_3N, int Np_3N, double *p_3N){

   
    // ---------------------------------------------- READ 3N --------------------------------------------

    bool print_content = false;

    char filename[300];
    std::strcpy(filename, file_path.c_str());

    hid_t   file_id_3, dataset_id_3;
    herr_t  ret_3;

    printf("read file %s\n", filename);

    file_id_3 = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // determine Np, Nq and Nalpha
    int N_h5_3[1];

    // N_p
    dataset_id_3 = H5Dopen(file_id_3, "Np", H5P_DEFAULT);
    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5_3);
    check_h5_read_call(ret_3);
    int Np_3N_3 = N_h5_3[0];
    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    // N_q
    dataset_id_3 = H5Dopen(file_id_3, "Nq", H5P_DEFAULT);
    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5_3);
    check_h5_read_call(ret_3);
    int Nq_3N_3 = N_h5_3[0];
    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    // N_alpha
    dataset_id_3 = H5Dopen(file_id_3, "Nalpha", H5P_DEFAULT);
    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5_3);
    check_h5_read_call(ret_3);
    int Jj_dim_3 = N_h5_3[0];
    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    // p_mesh
    size_t p_dst_size_3 =  sizeof( pmesh_table );
    pmesh_table p_h5_3[Np_3N_3];

    size_t p_dst_offset_3[3] = { HOFFSET( pmesh_table, index_ptable ),
                                 HOFFSET( pmesh_table, mesh_ptable ),
                                 HOFFSET( pmesh_table, weight_ptable )
                               };

    size_t p_dst_sizes_3[3] = { sizeof( p_h5_3[0].index_ptable ),
                                sizeof( p_h5_3[0].mesh_ptable ),
                                sizeof( p_h5_3[0].weight_ptable )
                              };

    double p_3N_3[Np_3N_3];
    double wp_3N_3[Np_3N_3];

    if(print_content){ printf("\np mesh:\n"); }

    ret_3 = H5TBread_table( file_id_3, "p mesh", p_dst_size_3, p_dst_offset_3, p_dst_sizes_3, p_h5_3 );
    check_h5_read_table_call(ret_3);
    for (int i = 0; i <= Np_3N_3 - 1; i++)
    {
        p_3N_3[i] = p_h5_3[i].mesh_ptable;
        wp_3N_3[i] = p_h5_3[i].weight_ptable;
       if(print_content){  printf("%d %.3f %.3f\n", i, p_3N_3[i], wp_3N_3[i]); }
    }

    // q_mesh
    double q_3N_3[Nq_3N_3];
    double wq_3N_3[Nq_3N_3];
    qmesh_table q_h5_3[Nq_3N_3];

    size_t q_dst_size_3 =  sizeof( qmesh_table );

    size_t q_dst_offset_3[3] = { HOFFSET( qmesh_table, index_qtable ),
                                 HOFFSET( qmesh_table, mesh_qtable ),
                                 HOFFSET( qmesh_table, weight_qtable )
                               };

    size_t q_dst_sizes_3[3] = { sizeof( q_h5_3[0].index_qtable ),
                                sizeof( q_h5_3[0].mesh_qtable ),
                                sizeof( q_h5_3[0].weight_qtable )
                              };

    if(print_content){ printf("q mesh:\n"); }

    ret_3 = H5TBread_table( file_id_3, "q mesh", q_dst_size_3, q_dst_offset_3, q_dst_sizes_3, q_h5_3 );
    check_h5_read_table_call(ret_3);
    for (int i = 0; i <= Nq_3N - 1; i++)
    {
        q_3N_3[i] = q_h5_3[i].mesh_qtable;
        wq_3N_3[i] = q_h5_3[i].weight_qtable;
        if(print_content){ printf("%d %.3f %.3f\n", i, q_3N_3[i], wq_3N_3[i]); }
    }

    if(print_content){ printf("Check mesh systems...\n"); }
    if ((Nq_3N != Nq_3N_3) || (Np_3N != Np_3N_3))
    {
        raise_error("Inconsistent mesh systems!");
        //printf("Inconsistent mesh systems!\n");
        //return 1;
    }

    for (int i = 0; i <= Np_3N - 1; i++)
    {
        if (fabs(p_3N[i] - p_3N_3[i]) > 1e-10)
        {
            raise_error("Inconsistent mesh systems!");
            //printf("Inconsistent mesh systems!\n");
            //return 1;
        }
    }

    for (int i = 0; i <= Nq_3N - 1; i++)
    {
        if (fabs(q_3N[i] - q_3N_3[i]) > 1e-10)
        {
            raise_error("Inconsistent mesh systems!");
            //printf("Inconsistent mesh systems!\n");
            //return 1;
        }
    }

    // PW table
    pw_table pw_h5_3[Jj_dim_3];
    size_t pw_dst_size_3 = sizeof( pw_table );
    size_t pw_dst_offset_3[10] = { HOFFSET( pw_table, alpha_pwtable ),
                                   HOFFSET( pw_table, L_pwtable ),
                                   HOFFSET( pw_table, S_pwtable ),
                                   HOFFSET( pw_table, J_pwtable ),
                                   HOFFSET( pw_table, T_pwtable ),
                                   HOFFSET( pw_table, l_pwtable ),
                                   HOFFSET( pw_table, twoj_pwtable ),
                                   HOFFSET( pw_table, twoJtotal_pwtable ),
                                   HOFFSET( pw_table, PARtotal_pwtable ),
                                   HOFFSET( pw_table, twoTtotal_pwtable )
                                 };

    size_t pw_dst_sizes_3[10] = { sizeof( pw_h5_3[0].alpha_pwtable ),
                                  sizeof( pw_h5_3[0].L_pwtable ),
                                  sizeof( pw_h5_3[0].S_pwtable ),
                                  sizeof( pw_h5_3[0].J_pwtable ),
                                  sizeof( pw_h5_3[0].T_pwtable ),
                                  sizeof( pw_h5_3[0].l_pwtable ),
                                  sizeof( pw_h5_3[0].twoj_pwtable ),
                                  sizeof( pw_h5_3[0].twoJtotal_pwtable ),
                                  sizeof( pw_h5_3[0].PARtotal_pwtable ),
                                  sizeof( pw_h5_3[0].twoTtotal_pwtable )
                                };

    int L12_Jj_3[Jj_dim_3];
    int S12_Jj_3[Jj_dim_3];
    int J12_Jj_3[Jj_dim_3];
    int T12_Jj_3[Jj_dim_3];
    int l3_Jj_3[Jj_dim_3];
    int two_j3_Jj_3[Jj_dim_3];

    ret_3 = H5TBread_table( file_id_3, "pw channels", pw_dst_size_3, pw_dst_offset_3, pw_dst_sizes_3, pw_h5_3 );
    check_h5_read_table_call(ret_3);
    if(print_content){ printf("\nNalpha = %d\n", Jj_dim_3);
                       printf("  i   L   S   J   T   l  2*j\n"); }
    for (int i = 0; i <= Jj_dim_3 - 1; i++)
    {
        L12_Jj_3[i] = pw_h5_3[i].L_pwtable;
        S12_Jj_3[i] = pw_h5_3[i].S_pwtable;
        J12_Jj_3[i] = pw_h5_3[i].J_pwtable;
        T12_Jj_3[i] = pw_h5_3[i].T_pwtable;
        l3_Jj_3[i] = pw_h5_3[i].l_pwtable;
        two_j3_Jj_3[i] = pw_h5_3[i].twoj_pwtable;
        if(print_content){ printf("%3d %3d %3d %3d %3d %3d %3d\n", i, L12_Jj_3[i], S12_Jj_3[i], J12_Jj_3[i], T12_Jj_3[i], l3_Jj_3[i], two_j3_Jj_3[i]); }
    }

    // matrix elements
    //H5D_layout_t layout_3;
    //hid_t cparms_3;
    //hid_t filespace_id_3;
    //hsize_t filespace_dims_3[6], chunk_dims_3[6];
    //herr_t filespace_status_n_3;
    //int filespace_rank_3, rank_chunk_3;

    // check if dimensions agree with matrix element data
    //MKL_INT64 P123_dim = Np_3N_3 * Nq_3N_3 * Jj_dim_3;
    //MKL_INT64 P123_dim_sq = P123_dim * P123_dim;

    if(print_content){ printf("Read matrix elements...\n"); }

    dataset_id_3 = H5Dopen(file_id_3, "matrix elements", H5P_DEFAULT);

    // Get dataset rank and dimension.
    /*filespace_id_3 = H5Dget_space(dataset_id_3);    // Get filespace handle first
    filespace_rank_3 = H5Sget_simple_extent_ndims(filespace_id_3);
    filespace_status_n_3 = H5Sget_simple_extent_dims(filespace_id_3, filespace_dims_3, NULL);

    printf("dataset rank: %d\ndataset dimensions: {%lu,%lu,%lu,%lu,%lu,%lu}\n\n",
           filespace_rank,
           (unsigned long)(filespace_dims[0]),
           (unsigned long)(filespace_dims[1]),
           (unsigned long)(filespace_dims[2]),
           (unsigned long)(filespace_dims[3]),
           (unsigned long)(filespace_dims[4]),
           (unsigned long)(filespace_dims[5])

          );

    cparms_3 = H5Dget_create_plist (dataset_id_3);
    layout_3 = H5Pget_layout (cparms_3);

    if (H5D_CHUNKED == layout_3)
    {
        //Get chunking information: rank and dimensions
        rank_chunk_3 = H5Pget_chunk(cparms_3, 6, chunk_dims_3);
        printf("chunk rank: %d\nchunk dimensions: {%lu,%lu,%lu,%lu,%lu,%lu}\n\n",
               rank_chunk,
               (unsigned long)(chunk_dims_3[0]),
               (unsigned long)(chunk_dims_3[1]),
               (unsigned long)(chunk_dims_3[2]),
               (unsigned long)(chunk_dims_3[3]),
               (unsigned long)(chunk_dims_3[4]),
               (unsigned long)(chunk_dims_3[5])
              );
    }*/

    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123);
    check_h5_read_call(ret_3);

    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    ret_3 = H5Fclose(file_id_3);
    check_h5_close_call(ret_3);
}

void get_h5_P123_dimensions(std::string file_path, int& Nalpha, int& Np, int& Nq){
    // ---------------------------------------------- READ 3N --------------------------------------------

    char filename[300];
    std::strcpy(filename, file_path.c_str());

    hid_t   file_id_3, dataset_id_3;
    herr_t  ret_3;

    printf("read file %s\n", filename);

    file_id_3 = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // determine Np, Nq and Nalpha
    int N_h5_3[1];

    // N_p
    dataset_id_3 = H5Dopen(file_id_3, "Np", H5P_DEFAULT);
    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5_3);
    check_h5_read_call(ret_3);
    int Np_3N_3 = N_h5_3[0];
    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    // N_q
    dataset_id_3 = H5Dopen(file_id_3, "Nq", H5P_DEFAULT);
    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5_3);
    check_h5_read_call(ret_3);
    int Nq_3N_3 = N_h5_3[0];
    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    // N_alpha
    dataset_id_3 = H5Dopen(file_id_3, "Nalpha", H5P_DEFAULT);
    ret_3 = H5Dread (dataset_id_3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5_3);
    check_h5_read_call(ret_3);
    int Jj_dim_3 = N_h5_3[0];
    ret_3 = H5Dclose(dataset_id_3);
    check_h5_close_call(ret_3);

    Nalpha = Jj_dim_3;
    Np     = Np_3N_3;
    Nq     = Nq_3N_3;
}

void calculate_permutation_operator(double* P123_array,
                                    int Nq, double* q_array, double* wq_array,
                                    int Np, double* p_array, double* wp_array,
                                    int Nx, double* x_array, double* wx_array,
                                    int Nalpha,
                                    int* L_2N_array,
                                    int* S_2N_array,
                                    int* J_2N_array,
                                    int* T_2N_array,
                                    int* l_3N_array,
                                    int* two_j_3N_array,
                                    int two_J_3N, int two_T_3N, int parity_3N){
    
    int P123_dim = Np * Nq * Nalpha;
    int P123_dim_sq = P123_dim * P123_dim;
    
    /* Notation change */
    int two_T = two_T_3N;
    int two_J = two_J_3N;
    int PAR   = parity_3N;
    int Nx_Gtilde = Nx;
    int Jj_dim = Nalpha;
    int Np_3N = Np;
    int Nq_3N = Nq;

    int* L12_Jj    = L_2N_array;
    int* S12_Jj    = S_2N_array;
    int* J12_Jj    = J_2N_array;
    int* T12_Jj    = T_2N_array;
    int* l3_Jj     = l_3N_array;
    int* two_j3_Jj = two_j_3N_array;

    double* p_3N = p_array;
    double* q_3N = q_array;
    /* End of notation change */

    int J12_max;
    double pmax_3N;
    double qmax_3N;

    // for ptilde and qtilde integrals
    int Npq;

    // for angular integrals in F_interpolate
    int Nz;

    // determine optimized Lmax: Lmax = max(get_L)+max(get_l)
    int max_L12 = 0;
    int max_l3 = 0;
    int max_J12 = 0;

    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
        if (J12_Jj[alpha] > max_J12) max_J12 = J12_Jj[alpha];
        if (L12_Jj[alpha] > max_L12) max_L12 = L12_Jj[alpha];
        if (l3_Jj[alpha] > max_l3) max_l3 = l3_Jj[alpha];
    }

    // for F_local_matrix prestorage and F_interpolate
    int lmax = GSL_MAX_INT(max_l3, max_L12) + 3; // for C4 it is possible to couple l=lmax with THREE Y_{1}^{mu}
    int l_interpolate_max = l_interpolate_max = 2 * (lmax - 3) + 3;

    std::cout << "lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";

    int Lmax = max_L12 + max_l3;
    int kLegendremax = 2 * max_L12;

    int two_jmax_Clebsch = 2 * lmax;
    int jmax_Clebsch = lmax;
    int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??

    // prestore Clebsch Gordan coefficients
    std::cout << "prestore ClebschGordan...\n";
    double *ClebschGordan_data = new double[(two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (2 * two_jmax_Clebsch + 1)];

    #pragma omp parallel
    {
        #pragma omp for collapse(3)

        for (int two_j1 = 0; two_j1 <= two_jmax_Clebsch; two_j1++){
            for (int two_j2 = 0; two_j2 <= two_jmax_Clebsch; two_j2++){
                for (int two_j3 = 0; two_j3 <= 2 * two_jmax_Clebsch; two_j3++){
                    for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1++){
                        for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2++){
                            ClebschGordan_data[((two_j1 + 1) * (two_j1 + 1) + two_m1 - two_j1 - 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (2 * two_jmax_Clebsch + 1)
                                               + ((two_j2 + 1) * (two_j2 + 1) + two_m2 - two_j2 - 1) * (2 * two_jmax_Clebsch + 1)
                                               + two_j3]
                                = ClebschGordan(two_j1, two_j2, two_j3, two_m1, two_m2, two_m1 + two_m2);
                        }
                    }
                }
            }
        }
    }

    std::cout << "prestore ClebschGordan_int...\n";
    double *ClebschGordan_int_data = new double[(jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (2 * jmax_Clebsch + 1)];

    #pragma omp parallel
    {
        #pragma omp for collapse(3)

        for (int j1 = 0; j1 <= jmax_Clebsch; j1++){
            for (int j2 = 0; j2 <= jmax_Clebsch; j2++){
                for (int j3 = 0; j3 <= 2 * jmax_Clebsch; j3++){
                    for (int m1 = -j1; m1 <= j1; m1++){
                        for (int m2 = -j2; m2 <= j2; m2++){
                            ClebschGordan_int_data[((j1 + 1) * (j1 + 1) + m1 - j1 - 1) * (jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (2 * jmax_Clebsch + 1)
                                                   + ((j2 + 1) * (j2 + 1) + m2 - j2 - 1) * (2 * jmax_Clebsch + 1)
                                                   + j3]
                                = ClebschGordan(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m1 + 2 * m2);
                        }
                    }
                }
            }
        }
    }

    int shat_dim = 1000;
    double *shat_data = new double[shat_dim];
    for (int j = 0; j <= shat_dim - 1; j++)
    {
        shat_data[j] = sqrt(2 * j + 1);
    }

    MKL_INT64 Pdim = Np_3N * Nq_3N * Jj_dim;
    double *P123_store = new double[Pdim * Pdim];

    for (MKL_INT64 i = 0; i <= Pdim * Pdim - 1; i++){
        P123_store[i] = 0.0;
    }

    // for angular integration in Gtilde
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

    std::cout << "prestore SixJ...\n";
    double *SixJ_array = new double[(two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)];
    std::cout << "SixJ_test (>0?):" << (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) << "\n";

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

    std::cout << "calculate matrix elements of permutation operator P123\n";
    
    pmax_3N = 0.0;
    qmax_3N = 0.0;

    for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
        for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
            for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
                Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
            }
        }
    }

    std::cout << "Gtilde\n";
    long int fullsize = Np_3N*Nq_3N*Nx_Gtilde*Jj_dim*(Jj_dim - 1);
    long int counter = 0;
    int frac_n, frac_o=0;
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
                            
                            counter += 1;
                            frac_n = (100*counter)/fullsize;
                            if (frac_n>frac_o){std::cout << frac_n << "%" << std::endl; frac_o=frac_n;}
                        }
                    }
                }
            }
        }
    }

    calculate_Ptilde_no_spline (P123_store, Pdim, (MKL_INT64) Np_3N, p_3N, (MKL_INT64) Nq_3N, q_3N, (MKL_INT64) Nx_Gtilde, x_Gtilde, wx_Gtilde, (MKL_INT64) Jj_dim, pmax_3N, qmax_3N, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, (MKL_INT64) Lmax, (MKL_INT64) max_L12, (MKL_INT64) max_l3, (MKL_INT64) two_J, (MKL_INT64) two_T, SixJ_array, two_jmax_SixJ, Gtilde_store);
}

void calculate_antisymmetrization_operator(std::string file_path, int &Np, int &Nq, int& Nalpha, double** A123, double* q_array, double* p_array){
    
    int D123_dim = Np * Nq * Nalpha;
    int D123_dim_sq = D123_dim * D123_dim;

    *A123 = new double [D123_dim_sq];
    
    //read_P123_bin_data_file(*A123, D123_dim_sq);
    //read_P123_csv_data_file(*A123);
    read_P123_h5_data_file(file_path, *A123, Nq, q_array, Np, p_array);
    
    /* Add square P123-term */
    #pragma omp parallel for 
    for (int idx=0; idx<D123_dim_sq; idx++){
        (*A123)[idx] *= 2;
    }

    
    /* Add identity term */
    #pragma omp parallel for 
    for (int alpha_idx=0; alpha_idx<Nalpha; alpha_idx++){
        int A123_diag_idx;
        for (int q_idx=0; q_idx<Nq; q_idx++){
            for (int p_idx=0; p_idx<Np; p_idx++){
                            
                A123_diag_idx = (int) (alpha_idx*Nq*Np + q_idx*Np + p_idx)*D123_dim +
                                       alpha_idx*Nq*Np + q_idx*Np + p_idx;

                (*A123)[A123_diag_idx] += 1;
            }
        }
    }

    #pragma omp parallel for 
    for (int idx=0; idx<D123_dim_sq; idx++){
        (*A123)[idx] /= 6;
    }
}
