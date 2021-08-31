
#ifndef TYPE_DEFS_H
#define TYPE_DEFS_H

#include <complex>
#include <string>

typedef unsigned int uint;
//~ typedef float floatType;
typedef double floatType;
typedef std::complex<floatType> cfloatType;
typedef std::complex<double> cdouble;

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
	bool        parallel_run;
	std::string potential_model;
	std::string subfolder;
	std::string grid_type;
	std::string parameter_walk;
	std::string energy_input_file;
	std::string average;
    std::string output_folder;
	std::string P123_folder;
} run_params;

#endif // TYPE_DEFS_H

