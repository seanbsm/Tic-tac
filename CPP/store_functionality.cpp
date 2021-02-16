#include "store_functionality.h"

/* Structs for compact storage */
typedef struct pw_table{
    int alpha_pwtable;
    int L_pwtable;
    int S_pwtable;
    int J_pwtable;
    int T_pwtable;
    int l_pwtable;
    int twoj_pwtable;
    int twoJtotal_pwtable;
    int PARtotal_pwtable;
    int twoTtotal_pwtable;
} pw_table;

typedef struct alpha_table{
    int alpha_alphaprime_index;
    int alpha_index;
    int alphaprime_index;
} alpha_table;


typedef struct pmesh_table{
    int index_ptable;
    double mesh_ptable;
} pmesh_table;


typedef struct qmesh_table{
    int index_qtable;
    double mesh_qtable;
} qmesh_table;

typedef struct Psparse_table{
    int index_row_Ptable;
    int index_col_Ptable;
    double value_Ptable;
} Psparse_table;

using namespace std;

// chop and store matrix elements in single precision
void store_sparse_matrix_elements_P123_h5 (double** P123_sparse_ptr_val_array,
                                           int**    P123_sparse_ptr_idx_array,
                                           int*     P123_sparse_ptr_dim_array,
                                           int  Np_3N, double *p_3N, int Nq_3N, double *q_3N,
                                           int  N_chn_3N,
                                           int* chn_3N_idx_array,
                                           int  Jj_dim,
                                           int *L12_Jj,
                                           int *S12_Jj,
                                           int *J12_Jj,
                                           int *T12_Jj,
                                           int *l3_Jj,
                                           int *two_j3_Jj,
                                           int *two_J_Jj,
                                           int *two_T_Jj,
                                           int *P_3N_array,
                                           std::string filename_in){
    MKL_INT64 D123_dim = Np_3N * Nq_3N * Jj_dim;
    MKL_INT64 D123_dim_sq = D123_dim * D123_dim;

    bool print_content = true;
    
    if (print_content){
        printf(" - Chop matrix elements...\n");
    }


    /* Convert filename to char-object */
    char filename[300];
    std::strcpy(filename, filename_in.c_str());

    int n_out_file;
    char out_file[200];
    n_out_file = sprintf(out_file, filename);

    if (print_content){
        cout << " - Write to: " << out_file << "\n";
    }

    hid_t cparms, dataset_id;
    herr_t status;

    hsize_t chunk_size = 10;
    int* fill_data = NULL;
    int compress  = 0;

    //hsize_t start_file [5];  /* Start of hyperslab in file   */
    //hsize_t count_file [5];  /* Block count in file          */
    //hsize_t start_mem  [5];  /* Start of hyperslab in memory */
    //hsize_t count_mem  [5];  /* Block count in memory        */

    hid_t file_id, group_id, mem_id;

    file_id = H5Fcreate(out_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Number of mesh points */
    int     N_h5  [1];
    hsize_t dim_N [1];

    dim_N[0] = 1;

    /* Np */
    N_h5[0]     = Np_3N+1;
    group_id    = H5Screate_simple(1, dim_N, NULL);
    dataset_id  = H5Dcreate(file_id, "Np", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status      = H5Dclose(dataset_id);
    status      = H5Sclose(group_id);

    /* Nq */
    N_h5[0]     = Nq_3N+1;
    group_id    = H5Screate_simple(1, dim_N, NULL);
    dataset_id  = H5Dcreate(file_id, "Nq", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status      = H5Dclose(dataset_id);
    status      = H5Sclose(group_id);

    /* Nalpha */
    N_h5[0]     = Jj_dim;
    group_id    = H5Screate_simple(1, dim_N, NULL);
    dataset_id  = H5Dcreate(file_id, "Nalpha", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status      = H5Dclose(dataset_id);
    status      = H5Sclose(group_id);
    printf(" - test 1 \n");
    /* N_chn_3N */
    N_h5[0]     = N_chn_3N;
    group_id    = H5Screate_simple(1, dim_N, NULL);
    dataset_id  = H5Dcreate(file_id, "N_chn_3N", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status      = H5Dclose(dataset_id);
    status      = H5Sclose(group_id);
    printf(" - test 2 \n");
    /* p-momentum mesh */
    pmesh_table p_dst_buf[Np_3N+1];

    /* Calculate the size and the offsets of our struct members in memory */
    size_t p_dst_size =  sizeof( pmesh_table );

    size_t p_dst_offset[2] = { HOFFSET( pmesh_table, index_ptable ),
                               HOFFSET( pmesh_table, mesh_ptable )
                             };

    /* Assign values to data structure */
    pmesh_table p_data[Np_3N+1];
    for (int i = 0; i < Np_3N+1; i++){
        p_data[i].index_ptable = i;
        p_data[i].mesh_ptable = p_3N[i];
        //p_data[i].weight_ptable = wp_3N[i];
    }

    const char *p_field_names[2]  = { "index", "mesh point" };

    hid_t p_field_type[2];
    p_field_type[0] = H5T_NATIVE_INT;
    p_field_type[1] = H5T_NATIVE_DOUBLE;

    group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

    status = H5TBmake_table("Mesh info p", group_id, "p mesh", 2, Np_3N+1,
                            p_dst_size, p_field_names, p_dst_offset, p_field_type,
                            chunk_size, fill_data, compress, p_data);

    /* Close the group */
    status = H5Gclose(group_id);

    /* q-momentum mesh */
    qmesh_table q_dst_buf[Nq_3N+1];

    /* Calculate the size and the offsets of our struct members in memory */
    size_t q_dst_size =  sizeof( qmesh_table );

    size_t q_dst_offset[2] = { HOFFSET( qmesh_table, index_qtable ),
                               HOFFSET( qmesh_table, mesh_qtable )
                             };

    /* Assign values to data structure */
    qmesh_table q_data[Nq_3N+1];
    for (int i = 0; i < Nq_3N+1; i++){
        q_data[i].index_qtable = i;
        q_data[i].mesh_qtable = q_3N[i];
        //q_data[i].weight_qtable = wq_3N[i];
    }

    const char *q_field_names[2]  = { "index", "mesh point"};

    hid_t q_field_type[2];
    q_field_type[0] = H5T_NATIVE_INT;
    q_field_type[1] = H5T_NATIVE_DOUBLE;

    group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

    status = H5TBmake_table("Mesh info q", group_id, "q mesh", 2, Nq_3N+1,
                            q_dst_size, q_field_names, q_dst_offset, q_field_type,
                            chunk_size, fill_data, compress, q_data);

    /* Close the group. */
    status = H5Gclose(group_id);
    
    /* PW quantum numbers */
    //pw_table pw_dst_buf[Jj_dim];

    /* Calculate the size and the offsets of our struct members in memory */
    size_t pw_dst_size = sizeof( pw_table );

    size_t pw_dst_offset[10] = { HOFFSET( pw_table, alpha_pwtable ),
                                 HOFFSET( pw_table, L_pwtable ),
                                 HOFFSET( pw_table, S_pwtable ),
                                 HOFFSET( pw_table, J_pwtable ),
                                 HOFFSET( pw_table, T_pwtable ),
                                 HOFFSET( pw_table, l_pwtable ),
                                 HOFFSET( pw_table, twoj_pwtable ),
                                 HOFFSET( pw_table, twoJtotal_pwtable ),
                                 HOFFSET( pw_table, PARtotal_pwtable ),
                                 HOFFSET( pw_table, twoTtotal_pwtable ),
                               };
                               
    /* Assign values to data structure */
    pw_table pw_data[Jj_dim];
    for (int i = 0; i <= Jj_dim - 1; i++){
        pw_data[i].alpha_pwtable = i;
        pw_data[i].L_pwtable = L12_Jj[i];
        pw_data[i].S_pwtable = S12_Jj[i];
        pw_data[i].J_pwtable = J12_Jj[i];
        pw_data[i].T_pwtable = T12_Jj[i];
        pw_data[i].l_pwtable = l3_Jj[i];
        pw_data[i].twoj_pwtable = two_j3_Jj[i];
        pw_data[i].twoJtotal_pwtable = two_J_Jj[i];
        pw_data[i].PARtotal_pwtable = P_3N_array[i];
        pw_data[i].twoTtotal_pwtable = two_T_Jj[i];
    }
    
    /* Define field information */
    const char *pw_field_names[10]  = { "index", "L_12", "S_12", "J_12", "T_12", "l_3", "2*j_3", "2*J_total", "PAR_total", "2*T_total" };

    hid_t pw_field_type[10];
    pw_field_type[0] = H5T_NATIVE_INT;
    pw_field_type[1] = H5T_NATIVE_INT;
    pw_field_type[2] = H5T_NATIVE_INT;
    pw_field_type[3] = H5T_NATIVE_INT;
    pw_field_type[4] = H5T_NATIVE_INT;
    pw_field_type[5] = H5T_NATIVE_INT;
    pw_field_type[6] = H5T_NATIVE_INT;
    pw_field_type[7] = H5T_NATIVE_INT;
    pw_field_type[8] = H5T_NATIVE_INT;
    pw_field_type[9] = H5T_NATIVE_INT;
    
    group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

    status = H5TBmake_table("Partial wave info", group_id, "pw channels", 10, Jj_dim,
                            pw_dst_size, pw_field_names, pw_dst_offset, pw_field_type,
                            chunk_size, fill_data, compress, pw_data);

    /* Close the group */
    status = H5Gclose(group_id);
    
    /* P123 block dimensions */
    hsize_t dims_file[1], dim_mem[1];
    dims_file[0] = N_chn_3N;

    dim_mem[0] = dims_file[0];

    group_id = H5Screate_simple(1, dims_file, NULL);

    hsize_t chunk_dim[] = {min(4, N_chn_3N)};
    cparms = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_chunk ( cparms, 1, chunk_dim);

    dataset_id = H5Dcreate(file_id, "P123 block dimensions", H5T_NATIVE_FLOAT, group_id, H5P_DEFAULT, cparms, H5P_DEFAULT);

    mem_id = H5Screate_simple(1, dim_mem, NULL);

    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_ptr_dim_array);

    status = H5Dclose(dataset_id);
    status = H5Sclose(mem_id);
    /* Close the group */
    status = H5Sclose(group_id);
    
     /* 3N channel indices */
    //hsize_t dims_file[1], dim_mem[1];
    dims_file[0] = N_chn_3N+1;

    dim_mem[0] = dims_file[0];

    group_id = H5Screate_simple(1, dims_file, NULL);

    //hsize_t chunk_dim[] = {min(4, N_chn_3N)};
    chunk_dim[0] = {min(4, N_chn_3N+1)};
    cparms = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_chunk ( cparms, 1, chunk_dim);

    dataset_id = H5Dcreate(file_id, "3N channel indices", H5T_NATIVE_FLOAT, group_id, H5P_DEFAULT, cparms, H5P_DEFAULT);

    mem_id = H5Screate_simple(1, dim_mem, NULL);

    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, chn_3N_idx_array);

    status = H5Dclose(dataset_id);
    status = H5Sclose(mem_id);
    /* Close the group */
    status = H5Sclose(group_id);

    /* Sparse matrix elements */
    for (int chn_3N_idx=0; chn_3N_idx<N_chn_3N; chn_3N_idx++){

        /* Sparse P123-matrix properties */
        double* P123_sparse_val_array = P123_sparse_ptr_val_array[chn_3N_idx];
        int*    P123_sparse_idx_array = P123_sparse_ptr_idx_array[chn_3N_idx];
        int     P123_sparse_dim       = P123_sparse_ptr_dim_array[chn_3N_idx];

        /* Dense P123-matrix decleration */
        //double* P123_dense_array_double_prec = NULL;
        //int     P123_dense_dim = Nalpha_block * Np_3N * Nq_3N;

        /* Convert sparse array to a dense format */
        //square_sparse_COO_to_dense_format_converter(P123_dense_dim,
        //                                            &P123_dense_array_double_prec,
        //                                            P123_sparse_val_array,
        //                                            P123_sparse_idx_array,
        //                                            P123_sparse_dim);

        /* Convert from double to float */
        //float* P123_dense_array = new float [P123_dense_dim];
        //for (int i=0; i<P123_dense_dim; i++){
        //    P123_dense_array[i] = P123_dense_array_double_prec[i];
        //}

        /* Calculate the size and the offsets of our struct members in memory */
        size_t Psparse_dst_size = sizeof( Psparse_table );

        size_t Psparse_dst_offset[3] = { HOFFSET( Psparse_table, index_row_Ptable ),
                                         HOFFSET( Psparse_table, index_col_Ptable ),
                                         HOFFSET( Psparse_table, value_Ptable )
                                        };
                                   
        /* Assign values to data structure */
        Psparse_table Psparse_data[P123_sparse_dim];
        for (int i=0; i<P123_sparse_dim; i++){
            Psparse_data[i].index_row_Ptable = P123_sparse_idx_array[i];
            Psparse_data[i].index_col_Ptable = P123_sparse_idx_array[i + P123_sparse_dim];
            Psparse_data[i].value_Ptable     = P123_sparse_val_array[i];
        }

        /* Define field information */
        const char *Psparse_field_names[3]  = { "row idx", "col idx", "P123 value"};

        hid_t Psparse_field_type[3];
        Psparse_field_type[0] = H5T_NATIVE_INT;
        Psparse_field_type[1] = H5T_NATIVE_INT;
        Psparse_field_type[2] = H5T_NATIVE_DOUBLE;

        /* Channel name */
        std::string chn_name_string = "matrix elements for 3N channel " + std::to_string(chn_3N_idx);
        char chn_name_char [200];
        n_out_file = sprintf(chn_name_char, chn_name_string.c_str());

        group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

        status = H5TBmake_table("P123 sparse info", group_id, chn_name_char, 3, P123_sparse_dim,
                                Psparse_dst_size, Psparse_field_names, Psparse_dst_offset, Psparse_field_type,
                                chunk_size, fill_data, compress, Psparse_data);

        /* Close the group */
        status = H5Gclose(group_id);
    }

    //status = H5Dclose(dataset_id);
    //status = H5Sclose(mem_id);
    //status = H5Sclose(group_id);
    status = H5Fclose(file_id);
}


// chop and store matrix elements in single precision
void store_matrix_elements_P123_h5 (double *P123_store,
                                    int Np_3N, double *p_3N, int Nq_3N, double *q_3N,
                                    int Jj_dim,
                                    int *L12_Jj,
                                    int *S12_Jj,
                                    int *J12_Jj,
                                    int *T12_Jj,
                                    int *l3_Jj,
                                    int *two_j3_Jj,
                                    int *two_J_Jj,
                                    int *two_T_Jj,
                                    std::string filename_in){

    MKL_INT64 D123_dim = Np_3N * Nq_3N * Jj_dim;
    MKL_INT64 D123_dim_sq = D123_dim * D123_dim;

    float *P123_float;
    P123_float = new float[D123_dim_sq];

    bool print_content = true;

    if (print_content){
        printf("chop matrix elements...\n");
    }


    double eps = 1e-10;
    for (MKL_INT64 alpha = 0; alpha <= Jj_dim - 1; alpha++)
    {
        for (MKL_INT64 q_index = 0; q_index <= Nq_3N - 1; q_index++)
        {
            for (MKL_INT64 p_index = 0; p_index <= Np_3N - 1; p_index++)
            {
                for (MKL_INT64 alpha_prime = 0; alpha_prime <= Jj_dim - 1; alpha_prime++)
                {
                    for (MKL_INT64 q_prime_index = 0; q_prime_index <= Nq_3N - 1; q_prime_index++)
                    {
                        for (MKL_INT64 p_prime_index = 0; p_prime_index <= Np_3N - 1; p_prime_index++)
                        {
                            MKL_INT64 index = (MKL_INT64) (alpha * Nq_3N * Np_3N + q_index * Np_3N + p_index) * D123_dim +
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
    }

    /* convert filename to char-object */
    char filename[300];
    std::strcpy(filename, filename_in.c_str());

    int n_out_file;
    char out_file[200];
    n_out_file = sprintf(out_file, filename);

    if (print_content){
        cout << "write to: " << out_file << "\n";
    }

    hid_t cparms, dataset_id;
    herr_t status;

    hsize_t chunk_size = 10;
    int *fill_data = NULL;
    int compress  = 0;

    hsize_t start_file[5];  /* Start of hyperslab in file*/
    hsize_t count_file[5];  /* Block count in file*/
    hsize_t start_mem[5];  /* Start of hyperslab in memory*/
    hsize_t count_mem[5];  /* Block count in memory*/

    hid_t file_id, group_id, mem_id;

    file_id = H5Fcreate(out_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

// NUMBER OF MESH POINTS
    int N_h5[1];
    hsize_t dim_N[1];
    dim_N[0] = 1;

// Np
    N_h5[0] = Np_3N+1;

    group_id = H5Screate_simple(1, dim_N, NULL);
    dataset_id = H5Dcreate(file_id, "Np", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status = H5Dclose(dataset_id);
    status = H5Sclose(group_id);

// Nq
    N_h5[0] = Nq_3N+1;

    group_id = H5Screate_simple(1, dim_N, NULL);
    dataset_id = H5Dcreate(file_id, "Nq", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status = H5Dclose(dataset_id);
    status = H5Sclose(group_id);

// Nalpha
    N_h5[0] = Jj_dim;

    group_id = H5Screate_simple(1, dim_N, NULL);
    dataset_id = H5Dcreate(file_id, "Nalpha", H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
    status = H5Dclose(dataset_id);
    status = H5Sclose(group_id);

// PMESH
    pmesh_table p_dst_buf[Np_3N+1];

// Calculate the size and the offsets of our struct members in memory
    size_t p_dst_size =  sizeof( pmesh_table );

    size_t p_dst_offset[3] = { HOFFSET( pmesh_table, index_ptable ),
                               HOFFSET( pmesh_table, mesh_ptable )
                             };

// assign values to data structure
    pmesh_table p_data[Np_3N+1];
    for (int i = 0; i < Np_3N+1; i++)
    {
        p_data[i].index_ptable = i;
        p_data[i].mesh_ptable = p_3N[i];
        //p_data[i].weight_ptable = wp_3N[i];
    }

    const char *p_field_names[2]  = { "index", "mesh point" };

    hid_t p_field_type[2];
    p_field_type[0] = H5T_NATIVE_INT;
    p_field_type[1] = H5T_NATIVE_DOUBLE;

    group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

    status = H5TBmake_table("Mesh info p", group_id, "p mesh", 2, Np_3N+1,
                            p_dst_size, p_field_names, p_dst_offset, p_field_type,
                            chunk_size, fill_data, compress, p_data);

// Close the group.
    status = H5Gclose(group_id);

// QMESH
    qmesh_table q_dst_buf[Nq_3N+1];

// Calculate the size and the offsets of our struct members in memory
    size_t q_dst_size =  sizeof( qmesh_table );

    size_t q_dst_offset[2] = { HOFFSET( qmesh_table, index_qtable ),
                               HOFFSET( qmesh_table, mesh_qtable )
                             };

// assign values to data structure
    qmesh_table q_data[Nq_3N];
    for (int i = 0; i < Nq_3N+1; i++)
    {
        q_data[i].index_qtable = i;
        q_data[i].mesh_qtable = q_3N[i];
        //q_data[i].weight_qtable = wq_3N[i];
    }

    const char *q_field_names[3]  = { "index", "mesh point", "mesh weight" };

    hid_t q_field_type[2];
    q_field_type[0] = H5T_NATIVE_INT;
    q_field_type[1] = H5T_NATIVE_DOUBLE;

    group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

    status = H5TBmake_table("Mesh info q", group_id, "q mesh", 2, Nq_3N+1,
                            q_dst_size, q_field_names, q_dst_offset, q_field_type,
                            chunk_size, fill_data, compress, q_data);


// Close the group.
    status = H5Gclose(group_id);

// PW quantum numbers
    pw_table pw_dst_buf[Jj_dim];

// Calculate the size and the offsets of our struct members in memory
    size_t pw_dst_size = sizeof( pw_table );

    size_t pw_dst_offset[10] = { HOFFSET( pw_table, alpha_pwtable ),
                                 HOFFSET( pw_table, L_pwtable ),
                                 HOFFSET( pw_table, S_pwtable ),
                                 HOFFSET( pw_table, J_pwtable ),
                                 HOFFSET( pw_table, T_pwtable ),
                                 HOFFSET( pw_table, l_pwtable ),
                                 HOFFSET( pw_table, twoj_pwtable ),
                                 HOFFSET( pw_table, twoJtotal_pwtable ),
                                 HOFFSET( pw_table, PARtotal_pwtable ),
                                 HOFFSET( pw_table, twoTtotal_pwtable ),
                               };

// assign values to data structure
    pw_table pw_data[Jj_dim];
    for (int i = 0; i <= Jj_dim - 1; i++)
    {
        pw_data[i].alpha_pwtable = i;
        pw_data[i].L_pwtable = L12_Jj[i];
        pw_data[i].S_pwtable = S12_Jj[i];
        pw_data[i].J_pwtable = J12_Jj[i];
        pw_data[i].T_pwtable = T12_Jj[i];
        pw_data[i].l_pwtable = l3_Jj[i];
        pw_data[i].twoj_pwtable = two_j3_Jj[i];
        pw_data[i].twoJtotal_pwtable = two_J_Jj[i];
        pw_data[i].PARtotal_pwtable = 1 - 2*( (L12_Jj[i]+l3_Jj[i])%2 ); // PAR= (-1)^(L+l)
        pw_data[i].twoTtotal_pwtable = two_T_Jj[i];
    }

// Define field information
    const char *pw_field_names[10]  = { "index", "L_12", "S_12", "J_12", "T_12", "l_3", "2*j_3", "2*J_total", "PAR_total", "2*T_total" };

    hid_t pw_field_type[10];
    pw_field_type[0] = H5T_NATIVE_INT;
    pw_field_type[1] = H5T_NATIVE_INT;
    pw_field_type[2] = H5T_NATIVE_INT;
    pw_field_type[3] = H5T_NATIVE_INT;
    pw_field_type[4] = H5T_NATIVE_INT;
    pw_field_type[5] = H5T_NATIVE_INT;
    pw_field_type[6] = H5T_NATIVE_INT;
    pw_field_type[7] = H5T_NATIVE_INT;
    pw_field_type[8] = H5T_NATIVE_INT;
    pw_field_type[9] = H5T_NATIVE_INT;

    group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

    status = H5TBmake_table("Partial wave info", group_id, "pw channels", 10, Jj_dim,
                            pw_dst_size, pw_field_names, pw_dst_offset, pw_field_type,
                            chunk_size, fill_data, compress, pw_data);

// Close the group.
    status = H5Gclose(group_id);

// matrix elements
    hsize_t dims_file[6], dim_mem[6];
    dims_file[0] = Jj_dim;
    dims_file[1] = Nq_3N;
    dims_file[2] = Np_3N;
    dims_file[3] = Jj_dim;
    dims_file[4] = Nq_3N;
    dims_file[5] = Np_3N;

    dim_mem[0] = dims_file[0];
    dim_mem[1] = dims_file[1];
    dim_mem[2] = dims_file[2];
    dim_mem[3] = dims_file[3];
    dim_mem[4] = dims_file[4];
    dim_mem[5] = dims_file[5];

    group_id = H5Screate_simple(6, dims_file, NULL);

    hsize_t chunk_dim[] = {min(4, Jj_dim), min(4, Nq_3N), min(4, Np_3N), min(4, Jj_dim), min(4, Nq_3N), min(4, Np_3N)};
    cparms = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_chunk ( cparms, 6, chunk_dim);

// COMPRESSION OPTIONS, probably not needed
//unsigned szip_options_mask;
//unsigned szip_pixels_per_block;
//szip_options_mask=H5_SZIP_NN_OPTION_MASK;
//szip_pixels_per_block=32;
//status = H5Pset_deflate( cparms, 6);
// there are problems for szip if chunk size is not appropriate
//H5Pset_szip (cparms, szip_options_mask, szip_pixels_per_block);

//    dataset_id = H5Dcreate(file_id, "matrix elements", H5T_NATIVE_DOUBLE, group_id, H5P_DEFAULT, cparms, H5P_DEFAULT);
    dataset_id = H5Dcreate(file_id, "matrix elements", H5T_NATIVE_FLOAT, group_id, H5P_DEFAULT, cparms, H5P_DEFAULT);

    mem_id = H5Screate_simple(6, dim_mem, NULL);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_float);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_store);

    status = H5Dclose(dataset_id);

    status = H5Sclose(mem_id);

    status = H5Sclose(group_id);

    status = H5Fclose(file_id);
}