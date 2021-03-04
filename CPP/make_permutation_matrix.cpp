
#include "make_permutation_matrix.h"

void read_P123_h5_data_file(std::string file_path,
							double* P123,
							int Np_3N, double *p_3N_3, double *wp_3N_3,
							int Nq_3N, double *q_3N_3, double *wq_3N_3){

   
	// ---------------------------------------------- READ 3N --------------------------------------------

	bool print_content = false;

	char filename[300];
	strcpy(filename, file_path.c_str());

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

	//for (int i = 0; i <= Np_3N - 1; i++)
	//{
	//    if (fabs(p_3N[i] - p_3N_3[i]) > 1e-10)
	//    {
	//        raise_error("Inconsistent mesh systems!");
	//        //printf("Inconsistent mesh systems!\n");
	//        //return 1;
	//    }
	//}

	//for (int i = 0; i <= Nq_3N - 1; i++)
	//{
	//    if (fabs(q_3N[i] - q_3N_3[i]) > 1e-10)
	//    {
	//        raise_error("Inconsistent mesh systems!");
	//        //printf("Inconsistent mesh systems!\n");
	//        //return 1;
	//    }
	//}

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

/* Add this if 6j takes a while */
void precalculate_Wigner_6j_symbols(){
}

void calculate_permutation_matrix_for_3N_channel(double** P123_val_dense_array,
												 double** P123_val_sparse_array,
												 int** P123_row_array,
												 int** P123_col_array,
												 int& P123_dim,
												 bool use_dense_format,
												 int  Nq, double* q_array, double* wq_array, int Np_per_WP, int Np_WP, double *p_array_WP_bounds,
												 int  Np, double* p_array, double* wp_array, int Nq_per_WP, int Nq_WP, double *q_array_WP_bounds,
												 int  Nx, double* x_array, double* wx_array,
												 int  Nalpha,
												 int* L_2N_array,
												 int* S_2N_array,
												 int* J_2N_array,
												 int* T_2N_array,
												 int* L_1N_array,
												 int* two_J_1N_array,
												 int  two_J_3N,
												 int  two_T_3N){
	
	/* Notation change */
	int two_J = two_J_3N;
	int two_T = two_T_3N;
	//int PAR   = parity_3N;
	int Nx_Gtilde = Nx;
	int Jj_dim = Nalpha;
	int Np_3N = Np;
	int Nq_3N = Nq;

	int* L12_Jj    = L_2N_array;
	int* S12_Jj    = S_2N_array;
	int* J12_Jj    = J_2N_array;
	int* T12_Jj    = T_2N_array;
	int* l3_Jj     = L_1N_array;
	int* two_j3_Jj = two_J_1N_array;

	bool print_content = true;
	bool print_Gtilde_progress = false;

	double* p_3N = p_array;
	double* q_3N = q_array;
	/* End of notation change */

	/* START OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */
	/* This code calculates the geometric function Gtilde_{alpha,alpha'}(p',q',x) as an array */

	int J12_max;
	double pmax_3N;
	double qmax_3N;

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

	if (print_content){
		std::cout << "lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";
	}

	int Lmax = max_L12 + max_l3;
	int kLegendremax = 2 * max_L12;

	int two_jmax_Clebsch = 2 * lmax;
	int jmax_Clebsch = lmax;
	int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??
	
	// prestore Clebsch Gordan coefficients
	//if (print_content){
	//    std::cout << "prestore ClebschGordan...\n";
	//}
	//double *ClebschGordan_data = new double[(two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (2 * two_jmax_Clebsch + 1)];
//
	//#pragma omp parallel
	//{
	//    #pragma omp for collapse(3)
//
	//    for (int two_j1 = 0; two_j1 <= two_jmax_Clebsch; two_j1++){
	//        for (int two_j2 = 0; two_j2 <= two_jmax_Clebsch; two_j2++){
	//            for (int two_j3 = 0; two_j3 <= 2 * two_jmax_Clebsch; two_j3++){
	//                for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1++){
	//                    for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2++){
	//                        ClebschGordan_data[((two_j1 + 1) * (two_j1 + 1) + two_m1 - two_j1 - 1) * (two_jmax_Clebsch + 1) * (two_jmax_Clebsch + 1) * (2 * two_jmax_Clebsch + 1)
	//                                           + ((two_j2 + 1) * (two_j2 + 1) + two_m2 - two_j2 - 1) * (2 * two_jmax_Clebsch + 1)
	//                                           + two_j3]
	//                            = ClebschGordan(two_j1, two_j2, two_j3, two_m1, two_m2, two_m1 + two_m2);
	//                    }
	//                }
	//            }
	//        }
	//    }
	//}

	//if (print_content){
	//    std::cout << "prestore ClebschGordan_int...\n";
	//}
	//double *ClebschGordan_int_data = new double[(jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (2 * jmax_Clebsch + 1)];
//
	//#pragma omp parallel
	//{
	//    #pragma omp for collapse(3)
//
	//    for (int j1 = 0; j1 <= jmax_Clebsch; j1++){
	//        for (int j2 = 0; j2 <= jmax_Clebsch; j2++){
	//            for (int j3 = 0; j3 <= 2 * jmax_Clebsch; j3++){
	//                for (int m1 = -j1; m1 <= j1; m1++){
	//                    for (int m2 = -j2; m2 <= j2; m2++){
	//                        ClebschGordan_int_data[((j1 + 1) * (j1 + 1) + m1 - j1 - 1) * (jmax_Clebsch + 1) * (jmax_Clebsch + 1) * (2 * jmax_Clebsch + 1)
	//                                               + ((j2 + 1) * (j2 + 1) + m2 - j2 - 1) * (2 * jmax_Clebsch + 1)
	//                                               + j3]
	//                            = ClebschGordan(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m1 + 2 * m2);
	//                    }
	//                }
	//            }
	//        }
	//    }
	//}
	
	//int shat_dim = 1000;
	//double *shat_data = new double[shat_dim];
	//for (int j = 0; j <= shat_dim - 1; j++)
	//{
	//    shat_data[j] = sqrt(2 * j + 1);
	//}
	
	//MKL_INT64 Pdim = Np_3N * Nq_3N * Jj_dim;
	//double *P123_store = new double[Pdim * Pdim];
	//
	//for (MKL_INT64 i = 0; i <= Pdim * Pdim - 1; i++){
	//    P123_store[i] = 0.0;
	//}
	
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

	if (print_content){
		std::cout << "prestore SixJ...\n";
	}
	double *SixJ_array = new double[(two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)];
	if (print_content){
		std::cout << "SixJ_test (>0?):" << (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) << "\n";
	}

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

	if (print_content){
		std::cout << "calculate permutation matrix elements\n";
	}

	pmax_3N = 0.0;
	qmax_3N = 0.0;

	/* WARNING: This loop ASSUMES two_J and two_T are already conserved between alpha and alphaprime */
	for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
		//int two_J = two_J_3N_array[alpha];
		//int two_T = two_T_3N_array[alpha];
		for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
			for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
				Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
			}
		}
	}

	if (print_content){
		std::cout << "Started working on Gtilde\n";
	}
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
						//int two_J = two_J_3N_array[alpha];
						for (MKL_INT64 alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
							Gtilde_store[alpha * Jj_dim * Np_3N * Nq_3N * Nx_Gtilde + alphaprime * Np_3N * Nq_3N * Nx_Gtilde + p_index * Nq_3N * Nx_Gtilde + q_index * Nx_Gtilde + x_index]
								= Gtilde_new (p_3N[p_index], q_3N[q_index], x_Gtilde[x_index], alpha, alphaprime, Jj_dim, Lmax, L12_Jj, l3_Jj, Atilde_store, two_J);

							counter += 1;
			
							if (print_Gtilde_progress){
								frac_n = (100*counter)/fullsize;
								if (frac_n>frac_o){std::cout << frac_n << "%" << std::endl; frac_o=frac_n;}
							}
						}
					}
				}
			}
		}
	}
	if (print_content){
		std::cout << "Finished working on Gtilde\n";
	}

	/* END OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */

	int P123_dense_dim = Np_WP * Nq_WP * Nalpha;

	/* Preallocate array if we use dense format. Otherwise (i.e. sparse) start with
	 * some reasonable guess (usually less than a percent), and expand if required. */
	int dense_array_size   = P123_dense_dim * P123_dense_dim;
	int sparse_step_length = 0;

	if (use_dense_format){
		*P123_val_dense_array = new double [dense_array_size];
		P123_dim        = P123_dense_dim;
	}
	else{
		if (dense_array_size>1000){
			sparse_step_length = dense_array_size/1000;
		}
		else{
			sparse_step_length = dense_array_size;
		}

		/* The sparse dimension is determined by counting, see the loops below */
		P123_dim = 0;

		*P123_val_sparse_array = new double [sparse_step_length];
		*P123_row_array = new int    [sparse_step_length];
		*P123_col_array = new int    [sparse_step_length];
	}

	int P123_mat_idx      = 0;
	int P123_row_idx      = 0;
	int P123_col_idx      = 0;
	int current_array_dim = sparse_step_length;
	/* <X_i'j'^alpha'| - loops (rows of P123) */
	for (int alphap_idx = 0; alphap_idx < Nalpha; alphap_idx++){
		if (print_content){
			printf("   - Working on row state %d/%d \n", alphap_idx, Nalpha);
		}
		for (int qp_idx_WP = 0; qp_idx_WP < Nq_WP; qp_idx_WP++){
			for (int pp_idx_WP = 0; pp_idx_WP < Np_WP; pp_idx_WP++){

				P123_row_idx = alphap_idx*Nq_WP*Np_WP + qp_idx_WP*Np_WP +  pp_idx_WP;

				/* |X_ij^alpha> - loops (columns of P123) */
				for (int alpha_idx = 0; alpha_idx < Nalpha; alpha_idx++){
					for (int q_idx_WP = 0; q_idx_WP < Nq_WP; q_idx_WP++){
						for (int p_idx_WP = 0; p_idx_WP < Np_WP; p_idx_WP++){

							P123_col_idx = alpha_idx*Nq_WP*Np_WP + q_idx_WP*Np_WP + p_idx_WP;
							
							double P123_element = calculate_P123_element_in_WP_basis (  alpha_idx,  p_idx_WP,  q_idx_WP, 
																					   alphap_idx, pp_idx_WP, qp_idx_WP, 
																					   Np_per_WP, p_array, wp_array,
																					   Nq_per_WP, q_array, wq_array,
																					   Nx,    x_array, wx_array,
																					   Np_WP, p_array_WP_bounds,
																					   Nq_WP, q_array_WP_bounds,
																					   Nalpha,
																					   Gtilde_store );
																					   
							if (use_dense_format){
								int P123_mat_idx              = (int) P123_row_idx*P123_dense_dim + P123_col_idx;
								(*P123_val_dense_array)[P123_mat_idx] = P123_element;
							}
							else if (P123_element!=0){  // For a sparse matrix this should enacted very few times
								/* Append to sparse value and index arrays */
								(*P123_val_sparse_array)[P123_dim] = P123_element;
								(*P123_row_array)[P123_dim]        = P123_row_idx;
								(*P123_col_array)[P123_dim]        = P123_col_idx;

								/* Increment sparse dimension (num of non-zero elements) */
								P123_dim += 1;

								/* If the dimension goes over the array dimension we increase the array size
								 * via a copy-paste-type routine, and increment the current array dimension */
								if ( P123_dim>=current_array_dim ){  // This should occur a small amount of the time
									increase_sparse_array_size(P123_val_sparse_array, current_array_dim, sparse_step_length);
									increase_sparse_array_size(P123_row_array,        current_array_dim, sparse_step_length);
									increase_sparse_array_size(P123_col_array,        current_array_dim, sparse_step_length);

									/* Increment sparse-array dimension */
									current_array_dim += sparse_step_length;
								}
							}
							else{ // If the element is zero and a sparse format is in use, simply move on
								continue;
							}
						}
					}
				}
			}
		}
	}

	/* Contract arrays to minimal size (number of non-zero elements) */
	if (use_dense_format==false){
		reduce_sparse_array_size(P123_val_sparse_array, current_array_dim, P123_dim);
		reduce_sparse_array_size(P123_row_array,        current_array_dim, P123_dim);
		reduce_sparse_array_size(P123_col_array,        current_array_dim, P123_dim);
	}

	/* Delete all temporary arrays */
	delete [] gtilde_array;
	delete [] Gtilde_store;
	delete [] Atilde_store;
	delete [] SixJ_array;
}

void calculate_permutation_matrices_for_all_3N_channels(double** P123_sparse_val_array,
														int**    P123_sparse_row_array,
														int**    P123_sparse_col_array,
														int&     P123_sparse_dim_array,
														int  Nq, double* q_array, double* wq_array, int Np_per_WP, int Np_WP, double *p_array_WP_bounds,
														int  Np, double* p_array, double* wp_array, int Nq_per_WP, int Nq_WP, double *q_array_WP_bounds,
														int  Nx, double* x_array, double* wx_array,
														int  Nalpha,
														int* L_2N_array,
														int* S_2N_array,
														int* J_2N_array,
														int* T_2N_array,
														int* L_1N_array,
														int* two_J_1N_array,
														int  two_J_3N,
														int  two_T_3N){
	
	bool print_content = true;

	/* Either use a dense-block or allocate space dynamically
	 * as non-zero elements are found. The 2nd approach is a bit slower,
	 * but necessary due to the typical size of the dense-block being unhandleable. */
	bool use_dense_format = false;
	
	/* Dense and sparse array dimensions */
	int     P123_array_dense_dim  = Nalpha * Np_WP * Nq_WP;
	int     P123_array_dim        = 0;

	/* Temporary dense matrix to be filled in subroutine if use_dense_format=true */
	double* P123_val_dense_array = NULL;

	if (print_content){
		printf(" - Begin calculating 3N-channel permutation matrix \n");
	}
	calculate_permutation_matrix_for_3N_channel(&P123_val_dense_array, P123_sparse_val_array, P123_sparse_row_array, P123_sparse_col_array, P123_array_dim, use_dense_format,
							  					Nq_WP*Nq_per_WP, q_array, wq_array, Np_per_WP, Np_WP, p_array_WP_bounds,
							  					Np_WP*Np_per_WP, p_array, wp_array, Nq_per_WP, Nq_WP, q_array_WP_bounds,
							  					Nx, x_array, wx_array,
							  					Nalpha,
							  					L_2N_array,
							  					S_2N_array,
							  					J_2N_array,
							  					T_2N_array,
							  					L_1N_array,
							  					two_J_1N_array,
												two_J_3N,
												two_T_3N);

	if (print_content){
		printf("   - Calculation finished. \n");
	}

	/* Convert dense-storage to sparse COO-format */
	if (use_dense_format){

		///* Old code snippet to check max matrix element, handy at times */ 
		//double max_element_dense = 0;
		//for (int idx=0; idx<P123_array_dim*P123_array_dim; idx++){
		//    if (abs(P123_val_dense_array[idx])>max_element_dense){
		//        max_element_dense = abs(P123_val_dense_array[idx]);
		//    } 
		//}
		//std::cout << "Max element dense: " << max_element_dense << std::endl;

		if (print_content){
			printf("   - Used dense format for calculation. Converting to sparse format \n");
		}
			
		square_dense_to_sparse_COO_format_converter(P123_array_dense_dim,
													P123_val_dense_array,
													P123_sparse_val_array,
													P123_sparse_row_array,
													P123_sparse_col_array,
													P123_sparse_dim_array);
	}
	else{
		P123_sparse_dim_array = P123_array_dim;
	}
		
	if (print_content){
		double P_123_subarray_density = (double) P123_sparse_dim_array / (P123_array_dense_dim*P123_array_dense_dim);

		printf(" - Dense dimension:   %dx%d \n", P123_array_dense_dim, P123_array_dense_dim);
		printf(" - Non-zero elements: %d \n",    P123_sparse_dim_array);
		printf(" - Density:           %.4f %% \n",  100*P_123_subarray_density);
	}
		
	delete [] P123_val_dense_array;
}