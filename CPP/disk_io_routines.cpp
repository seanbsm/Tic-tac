#include "disk_io_routines.h"

template <typename T>
std::string to_string_with_precision_and_sign(const T a_value, const int n = 6){
	std::string sgn = "";
	if (a_value>=0){sgn="+";}
    std::ostringstream out;
    out.precision(n);
    out << sgn << std::scientific
			   << a_value;
    return out.str();
}

using namespace std;

void open_file(std::ofstream &file,
			   std::string file_path,
			   bool rewrite_file){
	if (rewrite_file){
		/* Overwrite file */
		file.open(file_path);
	}
	else{
		/* Append to file */
		file.open(file_path, std::ios_base::app);
	}

	if (!file) {
        raise_error("Unable to open file " + file_path);
    }
}

void open_readfile(std::ifstream &file,
			   	  std::string file_path){

	file.open(file_path);

	if (!file) {
        raise_error("Unable to open file " + file_path);
    }
}

void read_parameter_sample_list(run_params run_parameters, std::vector<double>& parameter_vector, int& num_model_params, int& num_param_sets){

	std::string filename = run_parameters.parameter_file;

	/* Define file with input */
	std::ifstream infile(filename);
	
	std::string line;
	std::string delimiter = " ";
	double 	    value;

	std::vector<std::string> str_vector;
	
	bool first_line = true;

	/* Loop through lines in file */
	while (std::getline(infile, line)){
		/* Solution copied from https://www.delftstack.com/howto/cpp/cpp-split-string-by-space/
		 * and modified on 02/03/2022 */
		size_t pos = 0;
    	while ((pos = line.find(delimiter)) != string::npos){
			value = std::stod(line.substr(0, pos));
    	    parameter_vector.push_back(value);
    	    line.erase(0, pos + delimiter.length());
    	}
		/* Last value will not have a right-space */
		value = std::stod(line);
    	parameter_vector.push_back(value);
		if (first_line){
			num_model_params = parameter_vector.size();
			first_line = false;
		}
	}
	
	num_param_sets = parameter_vector.size() / num_model_params;
}

void read_input_energies(double*& energy_array,
						 int&	  num_energies,
						 std::string file_path){

	/* Only accept .txt-files for energy-input */
	size_t str_size = file_path.size();
	if (file_path.substr(str_size-4) != ".txt"){
		raise_error("Input energy file must be .txt");
	}

	/* Open file*/
	std::ifstream input_file;
	open_readfile(input_file, file_path);

	/* Read file */
	std::vector<double> energies_vector;
	double entry = 0;
	while (input_file >> entry){
        energies_vector.push_back(entry);
    }

	/* Close file */
	input_file.close();

	/* Copy input to vector */
	num_energies = energies_vector.size();
	energy_array = new double [num_energies];
	std::copy(energies_vector.begin(), energies_vector.end(), energy_array);
}

void store_sep_complex_matrix(double* re_matrix_array,
							  double* im_matrix_array,
                         	  size_t num_rows,
						 	  size_t num_cols,
						 	  size_t steplength,
						 	  std::string file_path,
						 	  bool rewrite_file,
						 	  std::string array_seperator_text){
	/* Open file*/
	std::ofstream output_file;
	open_file(output_file, file_path, rewrite_file);
	
	/* Fixes formatting of stored numbers */
	output_file << std::fixed
				<< std::showpos
				<< std::scientific
				//~ << std::right
				//~ << std::setw(14)
				<< std::setprecision(17);
			
	output_file << "# " << array_seperator_text << "\n";
	
	/* Append array-values */
	for (size_t r=0; r<num_rows; r++){
		for (size_t c=0; c<num_cols; c++){
			/* Append vector element */
			output_file << re_matrix_array[r*steplength+c] << " " << im_matrix_array[r*steplength+c] << "\n";
		}
	}
	
	/* Close writing session */
	output_file << std::endl;
	
	/* Close files */
	output_file.close();
}

void store_complex_vector_with_comments(cdouble* matrix_array,
										std::string* comment_array,
										size_t num_rows,
										std::string file_path,
										bool rewrite_file,
										std::string array_seperator_text){
	/* Open file*/
	std::ofstream output_file;
	open_file(output_file, file_path, rewrite_file);
	
	/* Fixes formatting of stored numbers */
	output_file << std::fixed
				<< std::showpos
				<< std::scientific
				//~ << std::right
				//~ << std::setw(14)
				<< std::setprecision(17);
			
	output_file << "# " << array_seperator_text << "\n";
	
	/* Append array-values */
	for (size_t r=0; r<num_rows; r++){
		/* Append vector element */
		output_file << matrix_array[r].real() << " " << matrix_array[r].imag()  << " " << comment_array[r] << "\n";
	}
	
	/* Close writing session */
	output_file << std::endl;
	
	/* Close files */
	output_file.close();
}

void store_complex_matrix(cdouble* matrix_array,
                          size_t num_rows,
						  size_t num_cols,
						  size_t steplength,
						  std::string file_path,
						  bool rewrite_file,
						  std::string array_seperator_text){
	/* Open file*/
	std::ofstream output_file;
	open_file(output_file, file_path, rewrite_file);
	
	/* Fixes formatting of stored numbers */
	output_file << std::fixed
				<< std::showpos
				<< std::scientific
				//~ << std::right
				//~ << std::setw(14)
				<< std::setprecision(17);
			
	output_file << "# " << array_seperator_text << "\n";
	
	/* Append array-values */
	for (size_t r=0; r<num_rows; r++){
		for (size_t c=0; c<num_cols; c++){
			/* Append vector element */
			output_file << matrix_array[r*steplength+c].real() << " " << matrix_array[r*steplength+c].imag() << "\n";
		}
	}
	
	/* Close writing session */
	output_file << std::endl;
	
	/* Close files */
	output_file.close();
}

void store_run_parameters(run_params run_parameters){
	
	std::string file_path = run_parameters.output_folder + "/" + "run_parameters.txt";
	
	/* Open file */
	std::ofstream output_file;
	open_file(output_file, file_path, true);
	
	/* Append input */
	output_file << create_input_printout_string(run_parameters);
	
	/* Close file */
	output_file.close();
}

void store_q_WP_kinematics_txt(size_t Nq_WP,
							   double* q_WP_array,
							   double* Eq_WP_array,
							   double* Tlab_WP_array,
							   double* q_WP_midpoint_array,
							   double* Eq_WP_midpoint_array,
							   double* Tlab_WP_midpoint_array,
						       std::string filename){
	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);

	result_file << "# First table of numbers are bin boundaries. \n";
	result_file << "# Second table of numbers are bin mid-points for energy WPs.\n";
	result_file << "# This means Ecm in the 2nd table is the mid-point of the corresponding Ecm boundaries.\n";
	result_file << "# \n";
	result_file << "# BOUNDARIES: \n";
	result_file << "# \n";
	
	std::string header1 = "# Index";
	std::string header2 = "q [MeV]";
	std::string header3 = "Ecm [MeV]";
	std::string header4 = "Tlab [MeV]";
	
	/* Headers */
	result_file << header1;
	result_file << std::right << std::setw(10 + 15 - header2.length()) << " " << header2;
	result_file << std::right << std::setw(10 + 15 - header3.length()) << " " << header3;
	result_file << std::right << std::setw(10 + 15 - header4.length()) << " " << header4 << "\n";

	/* Append array-values */
	for (size_t i=0; i<Nq_WP+1; i++){
		/* Fixes formatting of stored numbers */
		result_file << std::showpos
					<< std::right
                	<< std::setw(header2.length());
		/* Append vector element */
        result_file << i;
		/* Fixes formatting of stored numbers */
		result_file << std::right
					<< std::scientific
					<< std::setprecision(8);
		result_file << std::setw(10) << " " << q_WP_array[i];
		result_file << std::setw(10) << " " << Eq_WP_array[i];
		result_file << std::setw(10) << " " << Tlab_WP_array[i];
		result_file << "\n";
	}

	result_file << "# \n";
	result_file << "# BIN MID-POINTS: \n";
	result_file << "# \n";

	/* Headers */
	result_file << header1;
	result_file << std::right << std::setw(10 + 15 - header2.length()) << " " << header2;
	result_file << std::right << std::setw(10 + 15 - header3.length()) << " " << header3;
	result_file << std::right << std::setw(10 + 15 - header4.length()) << " " << header4 << "\n";

	/* Append array-values */
	for (size_t i=0; i<Nq_WP; i++){
		/* Fixes formatting of stored numbers */
		result_file << std::showpos
					<< std::right
                	<< std::setw(header2.length());
		/* Append vector element */
        result_file << i;
		/* Fixes formatting of stored numbers */
		result_file << std::right
					<< std::scientific
					<< std::setprecision(8);
		result_file << std::setw(10) << " " << q_WP_midpoint_array[i];
		result_file << std::setw(10) << " " << Eq_WP_midpoint_array[i];
		result_file << std::setw(10) << " " << Tlab_WP_midpoint_array[i];
		result_file << "\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void store_q_WP_boundaries_csv(fwp_statespace fwp_states,
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
	for (size_t i=0; i<fwp_states.Nq_WP+1; i++){
		/* Append vector element */
        result_file << fwp_states.q_WP_array[i] << "\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void restructure_breakup_data(std::complex<double>*    			 U_array,
							  solution_configuration 			 solve_config,
							  channel_os_indexing 				 chn_os_indexing,
							  run_params 						 run_parameters,
							  fwp_statespace 					 fwp_states,
							  swp_statespace 					 swp_states,
							  pw_3N_statespace 					 pw_states,
							  std::vector<std::string>&			 headers,
							  std::vector<std::complex<double>>& data_mat,
							  int&								 tot_num_rows,
							  int&								 tot_num_cols){
	
	/* Make local pointers for solution-configuration */
	double*	T_lab_array = solve_config.T_lab_array;
	double*	E_com_array = solve_config.E_com_array;
	double*	q_com_array = solve_config.q_com_array;
	/* Make local pointers & variables for pw-statespace */
	int	 two_J			= pw_states.two_J_3N_array[0];
	int  P_3N 			= pw_states.P_3N_array[0];
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	int* two_T_3N_array = pw_states.two_T_3N_array;
	/* Make local pointers & variables for SWP-statespace */
	int	   Np_WP   = swp_states.Np_WP;
	int	   Nq_WP   = swp_states.Nq_WP;
	double E_bound = swp_states.E_bound;
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;

	/* Create table headers 
	 * Only one column per unique alpha'-alpha breakup partial-wave channel.
	 * Each row is a (p'q')-(q) combination for given partial-wave channel */
	headers = {"real", "Tlab [MeV]" 		,
			   "real", "Ecm [MeV]"  		,
			   "real", "qcm [MeV] (init.)"	,
			   "real", "qlo [MeV] (init.)"	,
			   "real", "qup [MeV] (init.)"	};

	/* Locate:
	 * - number of unique alpha'-alpha breakup channels <alpha',p',q'|U|alpha,p,q>
	 * - number of rows needed for given input energy/momentum */
	std::map<int, int> alpha_map;		// Keeps track of columns
	std::map<int, int> row_idx_map;	// Maps correct element to given row & column
	std::vector<int> num_rows_vec = {0};
	int alpha_counter = 0;
	for (size_t q_idx=0; q_idx<num_q_com; q_idx++){
		int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[q_idx];
		int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[q_idx+1];

		std::map<int, int> row_num_map;	// Keeps track of rows in given column

		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			size_t idx_alpha_row = deuteron_idx_array[idx_d_row];
			std::string l = std::to_string(L_1N_array[idx_alpha_row]);
			std::string j = std::to_string(two_J_1N_array[idx_alpha_row]);

			/* Locate unique channels */
			for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
				size_t idx_alpha_col = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
				//size_t idx_q_col	 = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
				//size_t idx_p_col	 = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

				std::string Lp = std::to_string(L_2N_array[idx_alpha_col]);
				std::string Sp = std::to_string(S_2N_array[idx_alpha_col]);
				std::string Jp = std::to_string(J_2N_array[idx_alpha_col]);
				std::string lp = std::to_string(L_1N_array[idx_alpha_col]);
				std::string jp = std::to_string(two_J_1N_array[idx_alpha_col]);

				std::string label = "<"+Lp+","+Sp+","+Jp+","+lp+","+jp+"|U|"+l+","+j+">";
				
				/* Add new column if it doesn't exist */
				int key_alpha   = idx_alpha_col*pw_states.Nalpha + idx_alpha_row;
				int key_row_idx = idx_BU_chn*pw_states.Nalpha    + idx_alpha_row;
				if (alpha_map.find(key_alpha)==alpha_map.end()){ // New channel/column found
					headers.push_back("comp"); headers.push_back(label);
					headers.push_back("real"); headers.push_back("qlo [MeV]");
					headers.push_back("real"); headers.push_back("qup [MeV]");
					headers.push_back("real"); headers.push_back("plo [MeV]");
					headers.push_back("real"); headers.push_back("pup [MeV]");

					alpha_map[key_alpha]   = alpha_counter;
					alpha_counter 		  += 1;
				}

				if (row_num_map.find(key_alpha)==row_num_map.end()){ // New channel/column found
					row_idx_map[key_row_idx] = num_rows_vec[q_idx];
					row_num_map[key_alpha]   = num_rows_vec[q_idx] + 1;
				}
				else{// Channel/column alread exists, save element mapping and increment row-count
					row_idx_map[key_row_idx]  = row_num_map[key_alpha];	// Sets current number of rows in column/channel as (p,q)-pair index
					row_num_map[key_alpha]   += 1;
				}
			}
		}
		/* Read through map of rows per column and extract longest row*/
		int num_rows_for_q = 0;
		for (auto itr=row_num_map.begin(); itr!=row_num_map.end(); itr++){
			if (num_rows_for_q<itr->second){
				num_rows_for_q = itr->second;
			}
		}
		num_rows_vec.push_back(num_rows_for_q);
	}

	/* Read through map of rows per column and extract longest row*/
	tot_num_rows = num_rows_vec.back();
	tot_num_cols = headers.size()/2;

	/* Allocate data_mat and set to 0 */
	data_mat.resize(tot_num_cols* tot_num_rows);
	for (int i=0; i<data_mat.size(); i++){
		data_mat[i] = {0,0};
	}

	/* Loop over on-shell q-bins */
	for (size_t q_idx=0; q_idx<num_q_com; q_idx++){
		size_t q_WP_idx = q_com_idx_array[q_idx];
		int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[q_idx];
		int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[q_idx+1];

		double Tlab = T_lab_array[q_idx];
		double Ecm  = E_com_array[q_idx];
		double qcm  = q_com_array[q_idx];
		double qlo  = fwp_states.q_WP_array[q_WP_idx];
		double qhi  = fwp_states.q_WP_array[q_WP_idx+1];

		for (int row=num_rows_vec[q_idx]; row<num_rows_vec[q_idx+1]; row++){
			data_mat[row*tot_num_cols + 0] = {Tlab, 0};
			data_mat[row*tot_num_cols + 1] = {Ecm,  0};
			data_mat[row*tot_num_cols + 2] = {qcm,  0};
			data_mat[row*tot_num_cols + 3] = {qlo,  0};
			data_mat[row*tot_num_cols + 4] = {qhi,  0};
		}

		/* Pointer to either p_SWP_unco_array or p_SWP_coup_array,
		 * which is determined by whether the channel is coupled or not */
		double* e_SWP_array_ptr = NULL;

		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			size_t idx_alpha_row = deuteron_idx_array[idx_d_row];
			for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
				size_t idx_alpha_col = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
				size_t idx_q_col	 = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
				size_t idx_p_col	 = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

				size_t U_idx = idx_d_row*chn_os_indexing.num_BU_chns + idx_BU_chn;
				
				
				int L = L_2N_array[idx_alpha_col];
				int S = S_2N_array[idx_alpha_col];
				int J = J_2N_array[idx_alpha_col];
				int T = T_2N_array[idx_alpha_col];

				int two_T_3N = two_T_3N_array[idx_alpha_col];

				/* Detemine if this is a coupled channel.
				* !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_channel = false;
				bool state_1S0 = (S==0 && J==0 && L==0);
				bool coupled_via_L_2N = (run_parameters.tensor_force && L!=J && J!=0);
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_channel  = true;
				}

				if (coupled_channel){
					int chn_2N_idx = unique_2N_idx(L, S, J, T, coupled_channel, run_parameters);
					if (coupled_via_L_2N){
						if (L<J){
							e_SWP_array_ptr = &swp_states.e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1)];
						}
						else{
							e_SWP_array_ptr = &swp_states.e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1) + Np_WP+1];
						}
					}
					else if (coupled_via_T_3N){
						if (two_T_3N==1){
							e_SWP_array_ptr = &swp_states.e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1)];
						}
						else{
							e_SWP_array_ptr = &swp_states.e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1) + Np_WP+1];
						}
					}
					else{
						raise_error("Unknown coupling encountered in resolvent-calculation!");
					}
				}
				else{
					int chn_2N_idx = unique_2N_idx(L, S, J, T, coupled_channel, run_parameters);
					e_SWP_array_ptr = &swp_states.e_SWP_unco_array[chn_2N_idx * (Np_WP+1)];
				}

				//double qcm  = q_com_array[q_idx];
				double qlo  = fwp_states.q_WP_array[idx_q_col];
				double qup  = fwp_states.q_WP_array[idx_q_col+1];
				//double pcm  = s_com_array[q_idx];
				double plo  = com_energy_to_com_p_momentum(e_SWP_array_ptr[idx_p_col]);
				double pup  = com_energy_to_com_p_momentum(e_SWP_array_ptr[idx_p_col+1]);

				int key_alpha   = idx_alpha_col*pw_states.Nalpha + idx_alpha_row;
				int key_row_idx = idx_BU_chn*pw_states.Nalpha    + idx_alpha_row;

				int c = 5*alpha_map[key_alpha] + 5;
				int r = row_idx_map[key_row_idx];

				data_mat[r*tot_num_cols + c + 0] = U_array[U_idx];
				data_mat[r*tot_num_cols + c + 1] = {qlo, 0};
				data_mat[r*tot_num_cols + c + 2] = {qup, 0};
				data_mat[r*tot_num_cols + c + 3] = {plo, 0};
				data_mat[r*tot_num_cols + c + 4] = {pup, 0};
			}
		}
	}
}

void store_U_BU_matrix_elements_txt(std::complex<double>*    U_array,
								    solution_configuration solve_config,
								    channel_os_indexing chn_os_indexing,
								    run_params run_parameters,
									fwp_statespace fwp_states,
								    swp_statespace swp_states,
							        pw_3N_statespace pw_states,
							        std::string filename){
	
	/* Make local variable for run-parameters */
	std::string potential_model = run_parameters.potential_model;
	/* Make local pointers for solution-configuration */
	double*	T_lab_array = solve_config.T_lab_array;
	double*	E_com_array = solve_config.E_com_array;
	double*	q_com_array = solve_config.q_com_array;
	/* Make local pointers & variables for pw-statespace */
	int	 two_J			= pw_states.two_J_3N_array[0];
	int  P_3N 			= pw_states.P_3N_array[0];
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	int* two_T_3N_array = pw_states.two_T_3N_array;
	/* Make local pointers & variables for SWP-statespace */
	int	   Np_WP   = swp_states.Np_WP;
	int	   Nq_WP   = swp_states.Nq_WP;
	double E_bound = swp_states.E_bound;
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;

	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);

	std::string parity_sgn = "";
	if      (P_3N==+1){parity_sgn = "+";}
	else if	(P_3N==-1){parity_sgn = "-";}
	else			  {raise_error("Recieved illegal parity in store_U_matrix_elements_txt().");}

	result_file << std::setprecision(8);

	result_file << "# Breakup Nd-scattering U-matrix elements in a wave-packet, Jj-scheme representation. \n";
	result_file << "# \n";
	result_file << "# The calculations below were done for: \n";
	result_file << "# JP:          " << two_J << "/2" << parity_sgn << "\n";
	result_file << "# Potential:   " << potential_model << "\n";
	result_file << "# Np:          " << Np_WP << " \n";
	result_file << "# Nq:          " << Nq_WP << " \n";
	result_file << "# Particles:   " << "nd-scattering" << "\n";
	result_file << "# Deuteron BE: " << E_bound << " MeV \n";
	result_file << "# \n";
	result_file << "# Symbol definitions: \n";
	result_file << "# <L,S,J,2*l,2*j|U|2*l,2*j>:   U-matrix element in MeV with for channel <L S J l j|U|L=0 S=1 J=1 l j> (NOTE in wave-packet representation!)  \n";
	result_file << "# L:                           Pair orbital angular momentum  \n";
	result_file << "# S:                           Pair spin  \n";
	result_file << "# J:                           Pair total angular momentum  \n";
	result_file << "# l:                           Spectator orbital angular momentum  \n";
	result_file << "# j:                           Spectator total angular momentum  \n";
	result_file << "# Tlab:                        Laboratory scattering/kinetic energy  \n";
	result_file << "# Ecm:                         Centre-of-mass scattering/kinetic energy  \n";
	result_file << "# qcm:                         Spectator centre-of-mass momentum bin midpoint \n";
	result_file << "# qlo:                         Spectator centre-of-mass momentum bin lower boundary \n";
	result_file << "# qup:                         Spectator centre-of-mass momentum bin upper boundary \n";
	result_file << "# plo:                         Pair centre-of-mass momentum bin lower boundary \n";
	result_file << "# pup:                         Pair centre-of-mass momentum bin upper boundary \n";
	result_file << "# init.:                       Pre-scattering \n";
	result_file << "# \n";
	result_file << "# ################################################################################################### \n";
	//result_file << "#     Name    row-idx    col-idx     l'   2*j'      l    2*j  \n";

	result_file << std::right
				<< std::scientific
				<< std::setprecision(16);

	/* Widths of real- and complex-entry columns */
	int real_width = 24;
	int comp_width = 50;

	std::vector<std::string> headers;
	std::vector<std::complex<double>> 	 data_mat;
	int num_rows = 0;
	int num_cols = 0;
	
	restructure_breakup_data(U_array,
							 solve_config,
							 chn_os_indexing,
							 run_parameters,
							 fwp_states,
							 swp_states,
							 pw_states,
							 headers,
							 data_mat,
							 num_rows,
							 num_cols);
							 
	/* Write table headers */
	result_file << "#";
	for (size_t i=0; i<headers.size()/2; i++){
		if (headers[2*i]=="real"){
			result_file << std::right << std::setw(real_width) << headers[2*i+1];
		}
		else if (headers[2*i]=="comp"){
			result_file << std::right << std::setw(comp_width) << headers[2*i+1];
		}
		else{
			raise_error("Unknown argument in breakup amplitude disk i/o runtime!");
		}
	}
	result_file << "\n";
	
	/* Loop over data, row-wise */
	for (int r=0; r<num_rows; r++){
		result_file << " ";
		for (int c=0; c<num_cols; c++){
			if (headers[2*c]=="real"){
				result_file << std::right << std::setw(real_width) << data_mat[r*num_cols+c].real();
			}
			else if (headers[2*c]=="comp"){
				/* Append U-array element in Python numpy-style complex format */
				std::string U_string_real = to_string_with_precision_and_sign(data_mat[r*num_cols+c].real(), 16);
				std::string U_string_imag = to_string_with_precision_and_sign(data_mat[r*num_cols+c].imag(), 16);
				std::string U_string = U_string_real + U_string_imag + "j";

				result_file << std::right << std::setw(comp_width) << U_string;
			}
		}
		result_file << "\n";
	}
	result_file << "# ################################################################################################### \n";

	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void store_array(double* array, size_t array_length, std::string filename){
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
	for (size_t i=0; i<array_length; i++){
		/* Append vector element */
		//result_file << array[i] << "\n";
        result_file << array[i] << "\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}
void store_array(cdouble* array, size_t array_length, std::string filename){
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
	for (size_t i=0; i<array_length; i++){
		/* Append vector element */
		//result_file << array[i] << "\n";
        result_file << "(" << array[i].real() << array[i].imag() << "j)\n";
	}
	
	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void store_U_matrix_elements_txt(std::complex<double>*    U_array,
								 solution_configuration solve_config,
								 channel_os_indexing chn_os_indexing,
								 run_params run_parameters,
								 swp_statespace swp_states,
							     pw_3N_statespace pw_states,
							     std::string filename){
	
	/* Make local variable for run-parameters */
	std::string potential_model = run_parameters.potential_model;
	/* Make local pointers for solution-configuration */
	double*	T_lab_array = solve_config.T_lab_array;
	double*	E_com_array = solve_config.E_com_array;
	/* Make local pointers & variables for pw-statespace */
	int	 two_J			= pw_states.two_J_3N_array[0];
	int  P_3N 			= pw_states.P_3N_array[0];
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	/* Make local pointers & variables for SWP-statespace */
	int	   Np_WP   = swp_states.Np_WP;
	int	   Nq_WP   = swp_states.Nq_WP;
	double E_bound = swp_states.E_bound;
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;

	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);

	std::string parity_sgn = "";
	if      (P_3N==+1){parity_sgn = "+";}
	else if	(P_3N==-1){parity_sgn = "-";}
	else			  {raise_error("Recieved illegal parity in store_U_matrix_elements_txt().");}

	result_file << std::setprecision(8);

	result_file << "# Elastic Nd-scattering U-matrix elements in a wave-packet, Jj-scheme representation. \n";
	result_file << "# \n";
	result_file << "# The calculations below were done for: \n";
	result_file << "# JP:          " << two_J << "/2" << parity_sgn << "\n";
	result_file << "# Potential:   " << potential_model << "\n";
	result_file << "# Np:          " << Np_WP << " \n";
	result_file << "# Nq:          " << Nq_WP << " \n";
	result_file << "# Particles:   " << "nd-scattering" << "\n";
	result_file << "# Deuteron BE: " << E_bound << " MeV \n";
	result_file << "# \n";
	result_file << "# Symbol definitions: \n";
	result_file << "# Uij:   U-matrix element in MeV with row-idx i and col-idx j (NOTE in wave-packet representation!)  \n";
	result_file << "# l:     Spectator nucleon orbital angular momentum  \n";
	result_file << "# j:     Spectator nucleon total angular momentum  \n";
	result_file << "# Tlab:  Laboratory scattering/kinetic energy  \n";
	result_file << "# Ecm:   Centre-of-mass scattering/kinetic energy  \n";
	result_file << "# q_idx: Index of on-shell q-momentum bin correponding to Tlab/Ecm  \n";
	result_file << "# \n";
	result_file << "# ################################################################################################### \n";
	//result_file << "#     Name    row-idx    col-idx     l'   2*j'      l    2*j  \n";

	result_file << std::right
				<< std::scientific
				<< std::setprecision(16);

	/* Write headers for table 1 */
	std::vector<std::string> headers1 = {"Name", "row-idx", "col-idx", "l'", "2*j'", "l", "2*j"};
	result_file << "#";
	for (size_t idx_header=0; idx_header<headers1.size(); idx_header++){
		result_file << std::right << std::setw(10) << headers1[idx_header];
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

			std::string label = "U" + std::to_string(idx_d_row) + std::to_string(idx_d_col);

			result_file << " ";
			result_file << std::right << std::setw(10) << label;
			result_file << std::right << std::setw(10) << idx_d_row;
			result_file << std::right << std::setw(10) << idx_d_col;
			result_file << std::right << std::setw(10) << l_row;
			result_file << std::right << std::setw(10) << two_j_row;
			result_file << std::right << std::setw(10) << l_col;
			result_file << std::right << std::setw(10) << two_j_col;
			result_file << "\n";
		}
	}
	result_file << "# ################################################################################################### \n";

	/* Write headers for table 2 */
	std::vector<std::string> headers2 = {"Tlab [MeV]"  , "Ecm [MeV]"  , "q_idx"      ,
										 "U00 [MeV]"   , "U01 [MeV]"  , "U02 [MeV]"  ,
										 "U10 [MeV]"   , "U11 [MeV]"  , "U12 [MeV]"  ,
										 "U20 [MeV]"   , "U21 [MeV]"  , "U22 [MeV]"   };
	result_file << "#";
	for (size_t idx_header=0; idx_header<headers2.size(); idx_header++){
		if (idx_header<2){
			result_file << std::right << std::setw(24) << headers2[idx_header];
		}
		else if (idx_header==2){
			result_file << std::right << std::setw(10) << headers2[idx_header];
		}
		else{
			result_file << std::right << std::setw(50) << headers2[idx_header];
		}
	}
	result_file << "\n";

	/* Loop over on-shell q-bins */
	for (size_t q_idx=0; q_idx<num_q_com; q_idx++){
		size_t q_WP_idx = q_com_idx_array[q_idx];

		double T_lab = T_lab_array[q_idx];
		double E_com = E_com_array[q_idx];

		result_file << " ";
		result_file << std::right << std::setw(24) << T_lab;
		result_file << std::right << std::setw(24) << E_com;
		result_file << std::right << std::setw(10) << q_WP_idx;

		/* Loop over deuteron row-indices ("dp"=deuteron prime) */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			/* Loop over deuteron column-indices ("d"=deuteron) */
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				size_t U_idx = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + q_idx;
				
				std::string U_string_real = to_string_with_precision_and_sign(U_array[U_idx].real(), 16);
				std::string U_string_imag = to_string_with_precision_and_sign(U_array[U_idx].imag(), 16);
				std::string U_string = U_string_real + U_string_imag + "j";
				
				/* Append U-array element in Python numpy-style complex format */
    			result_file << std::right << std::setw(50) <<U_string;
			}
		}
		result_file << "\n";
	}

	result_file << "# ################################################################################################### \n";

	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}

void store_U_matrix_elements_csv(std::complex<double>*    U_array,
							     int* q_com_idx_array,	  size_t num_q_com,
					  		     int* deuteron_idx_array, size_t num_deuteron_states,
							     pw_3N_statespace pw_states,
							     std::string filename){
	
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;

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

void read_WP_boundaries_from_txt(double* WP_boundaries, int N_WP, std::string filename){
	std::ifstream infile(filename);
    std::string line;

    size_t idx_counter = 0;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        double val;
        if (!(iss >> val)){
			raise_error("Unknown format in input WP-boundary file:" + filename);
		}
        WP_boundaries[idx_counter] = val;
        idx_counter += 1;
    }
	if (idx_counter!=N_WP+1){
		raise_error("Program mismatch with number of wave-packets in input WP-boundaries file");
	}
}

void store_sparse_matrix_h5(double* sparse_val_array,
							int*    sparse_row_array,
							int*    sparse_col_array,
							size_t  sparse_dim,
							size_t  dense_dim,
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

	/* Write dimensions */
	if (print_content){
		printf("   - Writing matrix dimensions \n");
	}
	write_integer_to_h5((unsigned long long int) sparse_dim, "nnz", file_id);
	write_integer_to_h5((unsigned long long int) dense_dim,  "dim", file_id);
	
	/* Sparse matrix elements */
	if (print_content){
		printf("   - Writing sparse matrix elements and indices \n");
	}
	write_sparse_permutation_matrix_h5(sparse_val_array,
									   sparse_row_array,
									   sparse_col_array,
									   sparse_dim,
									   file_id);

	herr_t status = H5Fclose(file_id);
}

void read_sparse_matrix_h5(double** sparse_val_array,
						   int**    sparse_row_array,
						   int**    sparse_col_array,
						   size_t&  sparse_dim,
						   size_t&  dense_dim,
						   std::string filename_in,
						   bool    print_content){

	/* Convert filename_in to char-array */
	char filename[300];
	std::strcpy(filename, filename_in.c_str());

	if (print_content){
		printf(" - Read from: %s \n",  filename);
	}
	
	/* Read dimensions */
	unsigned long long int sparse_dim_temp = 0;
	read_ULL_integer_from_h5(sparse_dim_temp, "nnz", filename);
	sparse_dim = sparse_dim_temp;
	unsigned long long int dense_dim_temp = 0;
	read_ULL_integer_from_h5(dense_dim_temp, "dim", filename);
	dense_dim = dense_dim_temp;
	
	/* Read P123 sparse matrix elements and indices */
	*sparse_row_array = new int    [sparse_dim];
	*sparse_col_array = new int    [sparse_dim];
	*sparse_val_array = new double [sparse_dim];

	read_sparse_permutation_matrix_h5(*sparse_val_array,
									  *sparse_row_array,
									  *sparse_col_array,
									  sparse_dim,
									  filename);
}

void store_sparse_permutation_matrix_for_3N_channel_h5(double* P123_sparse_val_array,
													   int*    P123_sparse_row_array,
													   int*    P123_sparse_col_array,
													   size_t  P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   pw_3N_statespace pw_states,
													   std::string filename_in,
													   bool print_content){
	
	/* Make local pointers & variables */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	int  two_J_3N    	= pw_states.two_J_3N_array[0];
	int  P_3N   	    = pw_states.P_3N_array[0];

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
	write_ULL_integer_to_h5(P123_sparse_dim_temp, "P123_sparse_dim", file_id);

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
							  two_T_3N_array,
							  two_J_3N,
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
													   pw_3N_statespace pw_states,
													   std::string filename_in,
													   bool    print_content){
	
	/* Make local pointers & variables */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	int  two_J_3N    	= pw_states.two_J_3N_array[0];
	int  P_3N   	    = pw_states.P_3N_array[0];

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
	int two_T_3N_array_file [Nalpha];
	int two_J_3N_array_file [Nalpha];
	int P_3N_array_file     [Nalpha];
	read_PW_statespace_to_h5(Nalpha,
							 L_2N_array_file,
							 S_2N_array_file,
							 J_2N_array_file,
							 T_2N_array_file,
							 L_1N_array_file, 
							 two_J_1N_array_file,
							 two_T_3N_array_file,
							 two_J_3N_array_file,
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
			two_T_3N_array_file[i] != two_T_3N_array[i] ||
			two_J_3N_array_file[i] != two_J_3N ||
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

void write_integer_to_h5(int integer, const char* int_name, hid_t file_id){
	hid_t   group_id;
	hid_t   dataset_id;
	herr_t  status;

	int     N_h5  [1];
	hsize_t dim_N [1] = {1};

	/* Open file and create/write content correponding to variable-name int_name */
	N_h5[0]     = integer;
	group_id    = H5Screate_simple(1, dim_N, NULL);
	dataset_id  = H5Dcreate(file_id, int_name, H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5); check_h5_write_call(status);
	
	/* Close dataset and group */
	status      = H5Dclose(dataset_id); check_h5_close_call(status);
	status      = H5Sclose(group_id); check_h5_close_call(status);
}
void read_integer_from_h5(int& integer, const char* int_name, const char* filename){

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
	status = H5Dclose(dataset_id); check_h5_close_call(status);
	status = H5Fclose(file_id); check_h5_close_call(status);
}

void write_ULL_integer_to_h5(unsigned long long int integer, const char* int_name, hid_t file_id){
	hid_t   group_id;
	hid_t   dataset_id;
	herr_t  status;

	unsigned long long int N_h5  [1];
	hsize_t dim_N [1] = {1};

	/* Open file and create/write content correponding to variable-name int_name */
	N_h5[0]     = integer;
	group_id    = H5Screate_simple(1, dim_N, NULL);
	dataset_id  = H5Dcreate(file_id, int_name, H5T_NATIVE_ULLONG, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status      = H5Dwrite(dataset_id, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5); check_h5_write_call(status);
	
	/* Close dataset and group */
	status      = H5Dclose(dataset_id); check_h5_close_call(status);
	status      = H5Sclose(group_id); check_h5_close_call(status);
}
void read_ULL_integer_from_h5(unsigned long long int& integer, const char* int_name, const char* filename){

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
	unsigned long long int N_h5 [1];
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
	status = H5Dclose(dataset_id); check_h5_close_call(status);
	status = H5Fclose(file_id); check_h5_close_call(status);
}

void write_WP_boundaries_to_h5(double* WP_boundaries, int N_WP, const char* mesh_name, hid_t file_id){
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
void read_WP_boundaries_from_h5(double* WP_boundaries, int N_WP, const char* mesh_name, const char* filename){
	
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
							   int*  two_T_3N_array,
							   int   two_J_3N,
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
								 HOFFSET( pw_table, twoTtotal_pwtable ),
								 HOFFSET( pw_table, twoJtotal_pwtable ),
								 HOFFSET( pw_table, PARtotal_pwtable ),
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
		pw_data[i].twoTtotal_pwtable = two_T_3N_array[i];
		pw_data[i].twoJtotal_pwtable = two_J_3N;
		pw_data[i].PARtotal_pwtable  = P_3N;
	}
	
	/* Define field information */
	const char *pw_field_names[10]  = { "index", "L_12", "S_12", "J_12", "T_12", "l_3", "2*j_3", "2*T_total", "2*J_total", "PAR_total" };

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
							  int* two_T_3N_array,
							  int* two_J_3N_array,
							  int* P_3N_array,
							  const char* filename){
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
								 HOFFSET( pw_table, twoTtotal_pwtable ),
								 HOFFSET( pw_table, twoJtotal_pwtable ),
								 HOFFSET( pw_table, PARtotal_pwtable )
								};

	size_t pw_dst_sizes[10] = { sizeof( pw_h5[0].alpha_pwtable ),
								sizeof( pw_h5[0].L_pwtable ),
								sizeof( pw_h5[0].S_pwtable ),
								sizeof( pw_h5[0].J_pwtable ),
								sizeof( pw_h5[0].T_pwtable ),
								sizeof( pw_h5[0].l_pwtable ),
								sizeof( pw_h5[0].twoj_pwtable ),
								sizeof( pw_h5[0].twoTtotal_pwtable ),
								sizeof( pw_h5[0].twoJtotal_pwtable ),
								sizeof( pw_h5[0].PARtotal_pwtable )
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
		two_T_3N_array[i] = pw_h5[i].twoTtotal_pwtable;
		two_J_3N_array[i] = pw_h5[i].twoJtotal_pwtable;
		P_3N_array[i]     = pw_h5[i].PARtotal_pwtable;
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
									   const char*   filename){
	
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
