
#ifndef TYPE_DEFS_H
#define TYPE_DEFS_H

#include <complex>
#include <string>

typedef unsigned int uint;
//~ typedef float floatType;
typedef double floatType;
typedef std::complex<floatType> cfloatType;
typedef std::complex<double> cdouble;

typedef struct pw_3N_statespace{
	int  dim;
	int  J_2N_max;
	int* L_2N_array;
	int* S_2N_array;
	int* J_2N_array;
	int* T_2N_array;
	int* L_1N_array;
	int* two_J_1N_array;
	int* two_J_3N_array;
	int* two_T_3N_array;
	int* P_3N_array;
} pw_3N_statespace;

typedef struct run_params{
	int         two_J_3N_max;
	int         Np_WP;
	int         Nq_WP;
	int         J_2N_max;
	int         Nphi;
	int         Nx;
	double 		chebyshev_t;
	double 		chebyshev_s;
	int         Np_per_WP;
	int         Nq_per_WP;
	int			channel_idx;
	int			P123_omp_num_threads;
	int			max_TFC;
	bool        parallel_run;
	bool		P123_recovery;
	bool 		tensor_force;
	bool		isospin_breaking_1S0;
	bool 		mid_point_approximation;
	std::string potential_model;
	std::string subfolder;
	std::string p_grid_type;
	std::string p_grid_filename;
	std::string q_grid_type;
	std::string q_grid_filename;
	std::string parameter_walk;
	std::string energy_input_file;
	std::string average;
    std::string output_folder;
	std::string P123_folder;
} run_params;

#endif // TYPE_DEFS_H

