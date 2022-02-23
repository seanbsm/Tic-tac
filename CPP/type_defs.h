
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
	int  Nalpha;			// Number of partial waves, set in dynamical state-space construction (by construct_symmetric_pw_states)
	int  J_2N_max;			// Maximum pair-system total angular momentum
	int* L_2N_array;		// pair-nucleon state angular momentum
	int* S_2N_array;		// pair-nucleon state total spin
	int* J_2N_array;		// pair-nucleon state total angular momentum
	int* T_2N_array;		// pair-nucleon state total isospin
	int* L_1N_array;		// orbital-nucleon state angular momentum
	int* two_J_1N_array;	// orbital-nucleon state total angular momentum x2
	int* two_J_3N_array;	// three-nucleon state total angular momentum x2
	int* two_T_3N_array;	// three-nucleon state total isospin x2
	int* P_3N_array;		// three-nucleon state parity
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
	bool 		midpoint_approx;
	bool		calculate_and_store_P123;
	bool		solve_faddeev;
	bool		production_run;
	std::string potential_model;
	std::string subfolder;
	std::string p_grid_type;
	std::string p_grid_filename;
	std::string q_grid_type;
	std::string q_grid_filename;
	std::string parameter_walk;
	std::string energy_input_file;
    std::string output_folder;
	std::string P123_folder;
} run_params;

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

//typedef struct alpha_table{
//	int alpha_alphaprime_index;
//	int alpha_index;
//	int alphaprime_index;
//} alpha_table;
//
//
//typedef struct pmesh_table{
//	int index_ptable;
//	double mesh_ptable;
//} pmesh_table;
//
//
//typedef struct qmesh_table{
//	int index_qtable;
//	double mesh_qtable;
//} qmesh_table;
//
typedef struct bounds_mesh_table{
	int    index_table;
	double mesh_table;
} bounds_mesh_table;

//typedef struct Psparse_table{
//	int index_row_Ptable;
//	int index_col_Ptable;
//	double value_Ptable;
//} Psparse_table;

#endif // TYPE_DEFS_H

