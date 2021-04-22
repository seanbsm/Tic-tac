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

typedef struct bounds_mesh_table{
	int    index_table;
	double mesh_table;
} bounds_mesh_table;

typedef struct Psparse_table{
	int index_row_Ptable;
	int index_col_Ptable;
	double value_Ptable;
} Psparse_table;

using namespace std;

void store_q_WP_boundaries_csv(size_t Nq_WP, double* q_WP_array,
						   	   std::string filename){
	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);
	
	/* Fixes formatting of stored numbers */
	result_file << std::fixed
				<< std::showpos
				//~ << std::right
                //~ << std::setw(14)
				<< std::setprecision(8);
	
	/* Append array-values */
	for (size_t i=0; i<Nq_WP; i++){
		/* Append vector element */
        result_file << q_WP_array[i] << "\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void store_U_matrix_elements_csv(std::complex<double>* U_array,
							     int* q_com_idx_array,	  size_t num_q_com,
					  		     int* deuteron_idx_array, size_t num_deuteron_states,
							     int* L_1N_array, 
							     int* two_J_1N_array,
							     std::string filename){
	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);
	
	/* Fixes formatting of stored numbers */
	result_file << std::fixed
				<< std::showpos
				//~ << std::right
                //~ << std::setw(14)
				<< std::setprecision(8);
	
	result_file << "lp, two_jp, l, two_j";
	for (size_t q_idx=0; q_idx<num_q_com; q_idx++){
		size_t q_WP_idx = q_com_idx_array[q_idx];
		
		result_file << ", " << q_WP_idx;
	}
	result_file << "\n";

	/* Loop over deuteron row-indices ("dp"=deuteron prime) */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		size_t idx_alpha_row = deuteron_idx_array[idx_d_row];
		int l_row 	  = L_1N_array[idx_alpha_row];
		int two_j_row = two_J_1N_array[idx_alpha_row];

		/* Loop over deuteron column-indices ("d"=deuteron) */
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			size_t idx_alpha_col = deuteron_idx_array[idx_d_col];
			int l_col 	  = L_1N_array[idx_alpha_col];
			int two_j_col = two_J_1N_array[idx_alpha_col];

			result_file << l_row <<", "<< two_j_row <<", "<< l_col <<", "<< two_j_col;

			/* Loop over on-shell q-bins */
			for (size_t q_idx=0; q_idx<num_q_com; q_idx++){
				size_t q_WP_idx = q_com_idx_array[q_idx];
				
				size_t U_idx = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + q_idx;

				/* Append U-array element in Python numpy-style complex format */
        		result_file << ", (" << U_array[U_idx].real() << U_array[U_idx].imag() << "j)";
			}
			result_file << "\n";
		}
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void store_sparse_permutation_matrix_for_3N_channel_h5(double* P123_sparse_val_array,
													   int*    P123_sparse_row_array,
													   int*    P123_sparse_col_array,
													   size_t  P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   int     Nalpha,
													   int*    L_2N_array,
													   int*    S_2N_array,
													   int*    J_2N_array,
													   int*    T_2N_array,
													   int*    L_1N_array, 
													   int*    two_J_1N_array,
													   int     two_J_3N,
													   int     two_T_3N,
													   int     P_3N,
													   std::string filename_in,
													   bool print_content){

	if (print_content){
		printf("   - Setting up h5-file \n");
	}

	/* Convert filename_in to char-array */
	char filename[300];
	std::strcpy(filename, filename_in.c_str());

	if (print_content){
		cout << " - Write to: " << filename << "\n";
	}
	
	/* Create and open file */
	hid_t file_id = H5Fcreate(filename,
							  H5F_ACC_TRUNC,
							  H5P_DEFAULT,
							  H5P_DEFAULT);

	/* Write number of mesh points (Nalpha, Np_WP, Nq_WP, and P123_sparse_dim) */
	write_integer_to_h5(Nalpha,          "Nalpha",          file_id);
	write_integer_to_h5(Np_WP,           "Np_WP",           file_id);
	write_integer_to_h5(Nq_WP,           "Nq_WP",           file_id);

	unsigned long long int P123_sparse_dim_temp = P123_sparse_dim;
	write_integer_to_h5(P123_sparse_dim_temp, "P123_sparse_dim", file_id);

	/* Write p-momentum WP boundaries */
	if (print_content){
		printf("   - Writing p-momentum bins \n");
	}
	write_WP_boundaries_to_h5(p_WP_array, Np_WP, "p boundaries", file_id);
	/* Write q-momentum WP boundaries */
	if (print_content){
		printf("   - Writing q-momentum bins \n");
	}
	write_WP_boundaries_to_h5(q_WP_array, Nq_WP, "q boundaries", file_id);

	/* PW quantum numbers */
	if (print_content){
		printf("   - Writing partial-wave state space \n");
	}
	write_PW_statespace_to_h5(Nalpha,
							  L_2N_array,
							  S_2N_array,
							  J_2N_array,
							  T_2N_array,
							  L_1N_array, 
							  two_J_1N_array,
							  two_J_3N,
							  two_T_3N,
							  P_3N,
							  file_id);
	
	/* Sparse matrix elements */
	if (print_content){
		printf("   - Writing P123 sparse matrix elements and indices \n");
	}
	write_sparse_permutation_matrix_h5(P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim,
									   file_id);

	herr_t status = H5Fclose(file_id);
}

void read_sparse_permutation_matrix_for_3N_channel_h5(double** P123_sparse_val_array,
													   int**   P123_sparse_row_array,
													   int**   P123_sparse_col_array,
													   size_t& P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   int     Nalpha,
													   int*    L_2N_array,
													   int*    S_2N_array,
													   int*    J_2N_array,
													   int*    T_2N_array,
													   int*    L_1N_array, 
													   int*    two_J_1N_array,
													   int     two_J_3N,
													   int     two_T_3N,
													   int     P_3N,
													   std::string filename_in,
													   bool    print_content){

	/* Convert filename_in to char-array */
	char filename[300];
	std::strcpy(filename, filename_in.c_str());

	if (print_content){
		printf(" - Read from: %s \n",  filename);
	}
	
	/* Read number of mesh points (Nalpha, Np_WP, Nq_WP, and P123_sparse_dim) */
	int Nalpha_file = -1;
	int Np_WP_file  = -1;
	int Nq_WP_file  = -1;
	
	read_integer_from_h5(Nalpha_file,     "Nalpha",          filename);
	read_integer_from_h5(Np_WP_file,      "Np_WP",           filename);
	read_integer_from_h5(Nq_WP_file,      "Nq_WP",           filename);

	unsigned long long int P123_sparse_dim_temp = 0;
	read_ULL_integer_from_h5(P123_sparse_dim_temp, "P123_sparse_dim", filename);
	P123_sparse_dim = P123_sparse_dim_temp;

	/* Verify mesh points are equal to current program run, exit if not */
	if (Nalpha_file!=Nalpha || Np_WP_file!=Np_WP || Nq_WP_file!=Nq_WP){
		raise_error("File-read P123 state-space dimensions (Nalpha, Nq_WP, Np_WP) mismatch.");
	}

	/* Read p-momentum WP boundaries */
	double p_WP_array_file [Np_WP+1];
	read_WP_boundaries_from_h5(p_WP_array_file, Np_WP, "p boundaries", filename);
	/* Verify boundaries match current program run, exit if not */
	for (int i=0; i<Np_WP+1; i++){
		double p_boundary_prog = p_WP_array[i];
		double p_boundary_read = p_WP_array_file[i];
		if (p_boundary_read!=p_boundary_prog){
			raise_error("File-read P123 p-momentum boundaries mismatch.");
		}
	}

	/* Read q-momentum WP boundaries */
	double q_WP_array_file [Nq_WP+1];
	read_WP_boundaries_from_h5(q_WP_array_file, Nq_WP, "q boundaries", filename);
	/* Verify boundaries match current program run, exit if not */
	for (int i=0; i<Nq_WP+1; i++){
		double q_boundary_prog = q_WP_array[i];
		double q_boundary_read = q_WP_array_file[i];
		if (q_boundary_read!=q_boundary_prog){
			raise_error("File-read P123 q-momentum boundaries mismatch.");
		}
	}

	/* Read PW state space */
	int L_2N_array_file     [Nalpha];
	int S_2N_array_file     [Nalpha];
	int J_2N_array_file     [Nalpha];
	int T_2N_array_file     [Nalpha];
	int L_1N_array_file     [Nalpha];
	int two_J_1N_array_file [Nalpha];
	int two_J_3N_array_file [Nalpha];
	int two_T_3N_array_file [Nalpha];
	int P_3N_array_file     [Nalpha];
	read_PW_statespace_to_h5(Nalpha,
							 L_2N_array_file,
							 S_2N_array_file,
							 J_2N_array_file,
							 T_2N_array_file,
							 L_1N_array_file, 
							 two_J_1N_array_file,
							 two_J_3N_array_file,
							 two_T_3N_array_file,
							 P_3N_array_file,
							 filename);
	/* Verify PW statespace match current program run, exit if not */
	for (int i=0; i<Nalpha; i++){
		if (L_2N_array_file[i]     != L_2N_array[i] ||
			S_2N_array_file[i]     != S_2N_array[i] ||
			J_2N_array_file[i]     != J_2N_array[i] ||
			T_2N_array_file[i]     != T_2N_array[i] ||
			L_1N_array_file[i]     != L_1N_array[i] ||
			two_J_1N_array_file[i] != two_J_1N_array[i] ||
			two_J_3N_array_file[i] != two_J_3N ||
			two_T_3N_array_file[i] != two_T_3N ||
			P_3N_array_file[i]     != P_3N){
			raise_error("File-read P123 PW state-space mismatch.");
		}
	}

	/* Read P123 sparse matrix elements and indices */
	*P123_sparse_row_array = new int    [P123_sparse_dim];
	*P123_sparse_col_array = new int    [P123_sparse_dim];
	*P123_sparse_val_array = new double [P123_sparse_dim];

	read_sparse_permutation_matrix_h5(*P123_sparse_val_array,
									  *P123_sparse_row_array,
									  *P123_sparse_col_array,
									  P123_sparse_dim,
									  filename);
}

void write_integer_to_h5(int integer, char* int_name, hid_t file_id){
	hid_t   group_id;
	hid_t   dataset_id;
	herr_t  status;

	int     N_h5  [1];
	hsize_t dim_N [1] = {1};

	/* Open file and create/write content correponding to variable-name int_name */
	N_h5[0]     = integer;
	group_id    = H5Screate_simple(1, dim_N, NULL);
	dataset_id  = H5Dcreate(file_id, int_name, H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
	
	/* Close dataset and group */
	status      = H5Dclose(dataset_id);
	status      = H5Sclose(group_id);
}
void read_integer_from_h5(int& integer, char* int_name, char* filename){

	hid_t  file_id;
	hid_t  dataset_id;
	herr_t status;
	file_id = H5Fopen(filename,
					  H5F_ACC_RDONLY,
					  H5P_DEFAULT);

	/* Open file and find content correponding to variable-name int_name */
	dataset_id = H5Dopen(file_id,
						 int_name,
						 H5P_DEFAULT);

	/* Read from file into N_h5 */
	int N_h5 [1];
	status = H5Dread (dataset_id,
					  H5T_NATIVE_INT,
					  H5S_ALL,
					  H5S_ALL,
					  H5P_DEFAULT,
					  N_h5);
	check_h5_read_call(status);

	/* Write value to input integer */
	integer = N_h5[0];

	/* Close file */
	status = H5Dclose(dataset_id);
	check_h5_close_call(status);
	status = H5Fclose(file_id);
	check_h5_close_call(status);
}

void write_ULL_integer_to_h5(unsigned long long int integer, char* int_name, hid_t file_id){
	hid_t   group_id;
	hid_t   dataset_id;
	herr_t  status;

	unsigned long long int N_h5  [1];
	hsize_t dim_N [1] = {1};

	/* Open file and create/write content correponding to variable-name int_name */
	N_h5[0]     = integer;
	group_id    = H5Screate_simple(1, dim_N, NULL);
	dataset_id  = H5Dcreate(file_id, int_name, H5T_NATIVE_ULLONG, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status      = H5Dwrite(dataset_id, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5);
	
	/* Close dataset and group */
	status      = H5Dclose(dataset_id);
	status      = H5Sclose(group_id);
}
void read_ULL_integer_from_h5(unsigned long long& integer, char* int_name, char* filename){

	hid_t  file_id;
	hid_t  dataset_id;
	herr_t status;
	file_id = H5Fopen(filename,
					  H5F_ACC_RDONLY,
					  H5P_DEFAULT);

	/* Open file and find content correponding to variable-name int_name */
	dataset_id = H5Dopen(file_id,
						 int_name,
						 H5P_DEFAULT);

	/* Read from file into N_h5 */
	unsigned long long N_h5 [1];
	status = H5Dread (dataset_id,
					  H5T_NATIVE_ULLONG,
					  H5S_ALL,
					  H5S_ALL,
					  H5P_DEFAULT,
					  N_h5);
	check_h5_read_call(status);

	/* Write value to input integer */
	integer = N_h5[0];

	/* Close file */
	status = H5Dclose(dataset_id);
	check_h5_close_call(status);
	status = H5Fclose(file_id);
	check_h5_close_call(status);
}

void write_WP_boundaries_to_h5(double* WP_boundaries, int N_WP, char* mesh_name, hid_t file_id){
	/* Momentum mesh */
	bounds_mesh_table dst_buf[N_WP+1];

	/* Calculate the size and the offsets of our struct members in memory */
	size_t dst_size =  sizeof( bounds_mesh_table );

	size_t dst_offset[2] = { HOFFSET( bounds_mesh_table, index_table ),
							 HOFFSET( bounds_mesh_table, mesh_table )
						   };

	/* Assign values to data structure */
	bounds_mesh_table mesh_data[N_WP+1];
	for (int i = 0; i < N_WP+1; i++){
		mesh_data[i].index_table = i;
		mesh_data[i].mesh_table  = WP_boundaries[i];
	}

	const char *field_names[2]  = { "index", "mesh point" };

	hid_t field_type[2];
	field_type[0] = H5T_NATIVE_INT;
	field_type[1] = H5T_NATIVE_DOUBLE;

	hsize_t chunk_size = 10;
	int*    fill_data  = NULL;
	int     compress   = 0;

	hid_t   group_id;
	herr_t  status;

	group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

	char* mesh_info_name = NULL;
	if (mesh_name=="p boundaries"){
		mesh_info_name = {"Info p boundaries"};
	}
	else{
		mesh_info_name = {"Info q boundaries"};
	}

	status = H5TBmake_table(mesh_info_name, group_id, mesh_name, 2, N_WP+1,
							dst_size, field_names, dst_offset, field_type,
							chunk_size, fill_data, compress, mesh_data);

	/* Close the group */
	status = H5Gclose(group_id);
}
void read_WP_boundaries_from_h5(double* WP_boundaries, int N_WP, char* mesh_name, char* filename){
	
	hid_t  file_id;
	herr_t status;
	file_id = H5Fopen(filename,
					  H5F_ACC_RDONLY,
					  H5P_DEFAULT);
	
	/* p-mesh (fill p-momentum bin boundaries) */
	bounds_mesh_table p_h5[N_WP+1];

	size_t dst_size =  sizeof( bounds_mesh_table );
	
	size_t dst_offset[3] = { HOFFSET( bounds_mesh_table, index_table ),
							 HOFFSET( bounds_mesh_table, mesh_table )
							};

	size_t dst_sizes[3] = { sizeof( p_h5[0].index_table ),
							sizeof( p_h5[0].mesh_table )
						   };

	/* Read from file into p_h5 */
	status = H5TBread_table(file_id,
							mesh_name,
							dst_size,
							dst_offset,
							dst_sizes,
							p_h5);
	check_h5_read_table_call(status);

	/* Write from p_h5 into WP_boundaries */
	for (int i=0; i<N_WP+1; i++){
		WP_boundaries[i] = p_h5[i].mesh_table;
	}

	/* Close file */
	status = H5Fclose(file_id);
	check_h5_close_call(status);
}

void write_PW_statespace_to_h5(int   Nalpha,
							   int*  L_2N_array,
							   int*  S_2N_array,
							   int*  J_2N_array,
							   int*  T_2N_array,
							   int*  L_1N_array, 
							   int*  two_J_1N_array,
							   int   two_J_3N,
							   int   two_T_3N,
							   int   P_3N,
							   hid_t file_id){
	
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
	pw_table pw_data[Nalpha];
	for (int i=0; i<Nalpha; i++){
		pw_data[i].alpha_pwtable     = i;
		pw_data[i].L_pwtable         = L_2N_array[i];
		pw_data[i].S_pwtable         = S_2N_array[i];
		pw_data[i].J_pwtable         = J_2N_array[i];
		pw_data[i].T_pwtable         = T_2N_array[i];
		pw_data[i].l_pwtable         = L_1N_array[i];
		pw_data[i].twoj_pwtable      = two_J_1N_array[i];
		pw_data[i].twoJtotal_pwtable = two_J_3N;
		pw_data[i].PARtotal_pwtable  = P_3N;
		pw_data[i].twoTtotal_pwtable = two_T_3N;
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
	
	hsize_t chunk_size = 10;
	int*    fill_data  = NULL;
	int     compress   = 0;
	
	hid_t   group_id;
	herr_t  status;

	group_id = H5Gopen(file_id, "/", H5P_DEFAULT);

	status = H5TBmake_table("Partial wave info", group_id, "pw channels", 10, Nalpha,
							pw_dst_size, pw_field_names, pw_dst_offset, pw_field_type,
							chunk_size, fill_data, compress, pw_data);

	/* Close the group */
	status = H5Gclose(group_id);
}
void read_PW_statespace_to_h5(int  Nalpha,
							  int* L_2N_array,
							  int* S_2N_array,
							  int* J_2N_array,
							  int* T_2N_array,
							  int* L_1N_array, 
							  int* two_J_1N_array,
							  int* two_J_3N_array,
							  int* two_T_3N_array,
							  int* P_3N_array,
							  char* filename){
	hid_t  file_id;
	herr_t status;
	file_id = H5Fopen(filename,
					  H5F_ACC_RDONLY,
					  H5P_DEFAULT);
	
	/* Partial-wave table (fill PW arrays) */
	pw_table pw_h5[Nalpha];

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
								 HOFFSET( pw_table, twoTtotal_pwtable )
								};

	size_t pw_dst_sizes[10] = { sizeof( pw_h5[0].alpha_pwtable ),
								sizeof( pw_h5[0].L_pwtable ),
								sizeof( pw_h5[0].S_pwtable ),
								sizeof( pw_h5[0].J_pwtable ),
								sizeof( pw_h5[0].T_pwtable ),
								sizeof( pw_h5[0].l_pwtable ),
								sizeof( pw_h5[0].twoj_pwtable ),
								sizeof( pw_h5[0].twoJtotal_pwtable ),
								sizeof( pw_h5[0].PARtotal_pwtable ),
								sizeof( pw_h5[0].twoTtotal_pwtable )
							   };

	status = H5TBread_table(file_id,
							"pw channels",
							pw_dst_size,
							pw_dst_offset,
							pw_dst_sizes,
							pw_h5);
	check_h5_read_table_call(status);
	
	for (int i=0; i<Nalpha; i++){
		L_2N_array[i]     = pw_h5[i].L_pwtable;
		S_2N_array[i]     = pw_h5[i].S_pwtable;
		J_2N_array[i]     = pw_h5[i].J_pwtable;
		T_2N_array[i]     = pw_h5[i].T_pwtable;
		L_1N_array[i]     = pw_h5[i].l_pwtable;
		two_J_1N_array[i] = pw_h5[i].twoj_pwtable;
		two_J_3N_array[i] = pw_h5[i].twoJtotal_pwtable;
		P_3N_array[i]     = pw_h5[i].PARtotal_pwtable;
		two_T_3N_array[i] = pw_h5[i].twoTtotal_pwtable;
	}

	/* Close file */
	status = H5Fclose(file_id);
	check_h5_close_call(status);
}

void write_sparse_permutation_matrix_chunk_h5(double* P123_sparse_val_array,
											  int*    P123_sparse_row_array,
											  int*    P123_sparse_col_array,
											  int     P123_buffer_size,
											  hid_t   file_id){

}

void write_sparse_permutation_matrix_h5(double* P123_sparse_val_array,
										int*    P123_sparse_row_array,
										int*    P123_sparse_col_array,
										size_t  P123_sparse_dim,
										hid_t   file_id){

	hid_t       dataset_row;     /* dataset handle */
	hid_t       dataset_col;     /* dataset handle */
	hid_t       dataset_val;     /* dataset handle */
    hid_t       datatype_idx;	 /* handle */
	hid_t       datatype_val;	 /* handle */
	hid_t		dataspace;   	 /* handle */
    hsize_t     dimsf [1];       /* dataset dimensions */
	herr_t      status;
	
	/* Describe the size of the array and create the data space for fixed
     * size dataset */
	int RANK = 1;
    dimsf[0]  = P123_sparse_dim;
    dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
	/* Define datatype for the data in the file.
     * We will store little endian INT and DOUBLE numbers. */
    datatype_idx = H5Tcopy(H5T_NATIVE_INT);
	datatype_val = H5Tcopy(H5T_NATIVE_DOUBLE);
	
    status = H5Tset_order(datatype_idx, H5T_ORDER_LE);
	status = H5Tset_order(datatype_val, H5T_ORDER_LE);
	
	/* Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties. */
    dataset_row = H5Dcreate(file_id, "row indices",    datatype_idx, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	dataset_col = H5Dcreate(file_id, "column indices", datatype_idx, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	dataset_val = H5Dcreate(file_id, "values", 		   datatype_val, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	/* Write the data to the dataset using default transfer properties. */
    status = H5Dwrite(dataset_row, datatype_idx, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_row_array);
	status = H5Dwrite(dataset_col, datatype_idx, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_col_array);
	status = H5Dwrite(dataset_val, datatype_val, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_val_array);

	H5Sclose(dataspace);
    H5Tclose(datatype_idx);
    H5Dclose(dataset_row);
	H5Dclose(dataset_col);
	H5Dclose(dataset_val);

	///* Calculate the size and the offsets of our struct members in memory */
	//size_t Psparse_dst_size = sizeof( Psparse_table );
	//
	//size_t Psparse_dst_offset[3] = { HOFFSET( Psparse_table, index_row_Ptable ),
	//								 HOFFSET( Psparse_table, index_col_Ptable ),
	//								 HOFFSET( Psparse_table, value_Ptable )
	//								};
	//							   
	///* Assign values to data structure */
	//Psparse_table* Psparse_data = new Psparse_table [P123_sparse_dim];
	//for (int i=0; i<P123_sparse_dim; i++){
	//	Psparse_data[i].index_row_Ptable = P123_sparse_row_array[i];
	//	Psparse_data[i].index_col_Ptable = P123_sparse_col_array[i];
	//	Psparse_data[i].value_Ptable     = P123_sparse_val_array[i];
	//}
	//
	///* Define field information */
	//const char *Psparse_field_names[3]  = { "row idx", "col idx", "P123 value"};
	//
	//hid_t Psparse_field_type[3];
	//Psparse_field_type[0] = H5T_NATIVE_INT;
	//Psparse_field_type[1] = H5T_NATIVE_INT;
	//Psparse_field_type[2] = H5T_NATIVE_DOUBLE;
	//
	//hsize_t chunk_size = 0;
	//int max_chunk_size = 4*std::pow(2,30) / sizeof(double);
	//if (P123_sparse_dim<max_chunk_size){
	//	max_chunk_size = P123_sparse_dim;
	//}
	//else{
	//	max_chunk_size  = max_chunk_size / 3;
	//	max_chunk_size *= 3;
	//}
	//
	//chunk_size = max_chunk_size;
	//int*    fill_data  = NULL;
	//int     compress   = 0;
	//
	//hid_t   group_id;
	//herr_t  status;
	//
	//group_id = H5Gopen(file_id, "/", H5P_DEFAULT);
	//
	//status = H5TBmake_table("P123 sparse info", group_id, "sparse matrix", 3, P123_sparse_dim,
	//						Psparse_dst_size, Psparse_field_names, Psparse_dst_offset, Psparse_field_type,
	//						chunk_size, fill_data, compress, Psparse_data);
	//
	///* Close the group */
	//status = H5Gclose(group_id);
	//
	//delete [] Psparse_data;
}
void read_sparse_permutation_matrix_h5(double* P123_sparse_val_array,
									   int*    P123_sparse_row_array,
									   int*    P123_sparse_col_array,
									   size_t  P123_sparse_dim,
									   char*   filename){
	
	hid_t  file_id;
	hid_t  dataset_row;
	hid_t  dataset_col;
	hid_t  dataset_val;
	hid_t  datatype_idx;	 /* handle */
	hid_t  datatype_val;	 /* handle */
	herr_t status;

	/* Define datatype for the data in the file.
     * We will read little endian INT and DOUBLE numbers. */
    datatype_idx = H5Tcopy(H5T_NATIVE_INT);
	datatype_val = H5Tcopy(H5T_NATIVE_DOUBLE);

	/* open file */
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	/* Open file and find content correponding to variable-name int_name */
	dataset_row = H5Dopen(file_id, "row indices", 	 H5P_DEFAULT);
	dataset_col = H5Dopen(file_id, "column indices", H5P_DEFAULT);
	dataset_val = H5Dopen(file_id, "values", 		 H5P_DEFAULT);

	/* Read from file into N_h5 */
	status = H5Dread (dataset_row, datatype_idx, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_row_array);
	status = H5Dread (dataset_col, datatype_idx, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_col_array);
	status = H5Dread (dataset_val, datatype_val, H5S_ALL, H5S_ALL, H5P_DEFAULT, P123_sparse_val_array);

	/* Close file */
	status = H5Dclose(dataset_row); check_h5_close_call(status);
	status = H5Dclose(dataset_col); check_h5_close_call(status);
	status = H5Dclose(dataset_val); check_h5_close_call(status);
	status = H5Fclose(file_id);     check_h5_close_call(status);

	//hid_t  file_id;
	//herr_t status;
	//file_id = H5Fopen(filename,
	//				  H5F_ACC_RDONLY,
	//				  H5P_DEFAULT);
	//
	///* Calculate the size and the offsets of our struct members in memory */
	////Psparse_table P123_h5[P123_sparse_dim];
	//Psparse_table* P123_h5 = new Psparse_table [P123_sparse_dim];
	//
	//size_t P123_dst_size = sizeof( Psparse_table );
	//
	//size_t P123_dst_offset[3] = { HOFFSET( Psparse_table, index_row_Ptable ),
	//							  HOFFSET( Psparse_table, index_col_Ptable ),
	//							  HOFFSET( Psparse_table, value_Ptable )
	//							};
	//size_t P123_dst_sizes[3] =  { sizeof( P123_h5[0].index_row_Ptable ),
	//							  sizeof( P123_h5[0].index_col_Ptable ),
	//							  sizeof( P123_h5[0].value_Ptable )
	//							};
	//							
	///* Read from file into p_h5 */
	//status = H5TBread_table(file_id,
	//						"sparse matrix",
	//						P123_dst_size,
	//						P123_dst_offset,
	//						P123_dst_sizes,
	//						P123_h5);
	//check_h5_read_table_call(status);
	//
	///* Close file */
	//status = H5Fclose(file_id);
	//check_h5_close_call(status);
	//
	///* Write read data into argument-arrays */
	//for (int i=0; i<P123_sparse_dim; i++){
	//	P123_sparse_row_array[i] = P123_h5[i].index_row_Ptable;
	//	P123_sparse_col_array[i] = P123_h5[i].index_col_Ptable;
	//	P123_sparse_val_array[i] = P123_h5[i].value_Ptable;
	//}
	//
	//delete [] P123_h5;
}
