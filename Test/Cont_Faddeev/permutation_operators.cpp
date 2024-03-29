
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

void read_P123_h5_data_file(double* P123,
                            int Nq_3N, double *q_3N,
                            int Np_3N, double *p_3N,
                            int Nalpha_3N, int* L12_Jj_3, int* S12_Jj_3, int* J12_Jj_3, int* T12_Jj_3, int* l3_Jj_3, int* two_j3_Jj_3){

   
    // ---------------------------------------------- READ 3N --------------------------------------------

    bool print_content = false;

    char filename[300];

    hid_t   file_id_3, dataset_id_3;
    herr_t  ret_3;

   // double *P123_store = new double[V3N_Jj_dim_sq];

    //for (MKL_INT64 i = 0; i <= V3N_Jj_dim_sq - 1; i++) {
    //    P123[i] = 0.0;
    //}

    int i_file3 = sprintf(filename, "../../../Data/3N_permutation_operator/P123_files/P123_Jmax_1_Np_32_Nq_30.h5");
    if (i_file3 < 0){
        raise_error("Couldn't locate P123-matrix h5-file.");
    }

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

    //int L12_Jj_3[Jj_dim_3];
    //int S12_Jj_3[Jj_dim_3];
    //int J12_Jj_3[Jj_dim_3];
    //int T12_Jj_3[Jj_dim_3];
    //int l3_Jj_3[Jj_dim_3];
    //int two_j3_Jj_3[Jj_dim_3];

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

void read_P123_csv_data_file(double* P123){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    std::string P123_file_path = "../../../Data/3N_permutation_operator/Faddeev_permutation_operator_benchmark/P123_matrix_elements.csv";
    //std::string P123_file_path = "text_csv.csv";
    
    // Create an input filestream
    std::ifstream data_file(P123_file_path);

    // Make sure the file is open
    if(!data_file.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line;
    double val;

    int incrementing_idx = 0;

    // Read data, line by line
    while(std::getline(data_file, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each integer
        while(ss >> val){
            
            // Add the current integer to the 'colIdx' column's values vector
            P123[incrementing_idx] = val;
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
            // Increment the column index
            colIdx++;
        }
        incrementing_idx += 1;
    }

    // Close file
    data_file.close();
}

void read_P123_bin_data_file(double* P123, int P123_array_size){   

    std::string P123_file_path = "../../../Data/3N_permutation_operator/Faddeev_permutation_operator_benchmark/P123_matrix_elements_2.bin";

    std::ifstream infile(P123_file_path.c_str(), std::ios::in | std::ios::binary);
    infile.seekg(0, std::ios::end); 
    int size = infile.tellg();
    std::cout << size/sizeof(float) << " " << P123_array_size << std::endl;

    float binData [P123_array_size];
    //or float binData[numRows][numCols]; if you knew it can fit on the stack 
    //and your compiler has variable sized arrays

    //read into the array
    infile.read(reinterpret_cast<char *>(binData), sizeof(float) * P123_array_size);

    /*std::ifstream file("../../../Data/3N_permutation_operator/Faddeev_permutation_operator_benchmark/P123_matrix_elements.bin", std::ios::in | std::ios::binary | std::ios::ate);
    file.seekg(0, std::ios::end); 
    int size = file.tellg();  
    file.seekg(0, std::ios::beg); 
    char* memblock = new char[size];
    file.read(memblock, size);
    file.close();
    float* element_array = (float*) memblock; //reinterpret as float
    
    //double g = element_array[P123_array_size-1];
    std::cout << size << std::endl;*/

    //#pragma omp parallel for
    /*std::cout << "start reading P123" << std::endl;
    for (int idx=0; idx<P123_array_size; idx++){
        P123[idx] = (double)element_array[idx];//reinterpret as double
    }
    std::cout << "end reading P123" << std::endl;*/
    
    //delete [] memblock;
    //delete [] element_array;
}

void get_h5_P123_dimensions(int& Nalpha, int& Np, int& Nq){
    // ---------------------------------------------- READ 3N --------------------------------------------

    char filename[300];

    hid_t   file_id_3, dataset_id_3;
    herr_t  ret_3;

    int i_file3 = sprintf(filename, "../../../Data/3N_permutation_operator/P123_files/P123_Jmax_1_Np_32_Nq_30.h5");
    if (i_file3 < 0){
        raise_error("Couldn't locate P123-matrix h5-file.");
    }

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

void calculate_antisymmetrization_operator(int &Np, int &Nq, int& Nalpha, double** A123, double* q_array, double* p_array){
    
    int D123_dim = Np * Nq * Nalpha;
    int D123_dim_sq = D123_dim * D123_dim;

    *A123 = new double [D123_dim_sq];
    
    //read_P123_bin_data_file(*A123, D123_dim_sq);
    //read_P123_csv_data_file(*A123);
    //read_P123_h5_data_file(*A123, Nq, q_array, Np, p_array);
    
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
