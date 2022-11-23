
#ifndef TYPE_DEFS_H
#define TYPE_DEFS_H

#include <complex>
#include <string>

typedef unsigned int uint;
//~ typedef float floatType;
typedef double floatType;
typedef std::complex<floatType> cfloatType;
typedef std::complex<double> cdouble;
typedef unsigned long long int ull_int;

typedef struct pw_3N_statespace{
	int  Nalpha;				// Number of partial waves, set in dynamical state-space construction (by construct_symmetric_pw_states)
	int  J_2N_max;				// Maximum pair-system total angular momentum
	int* L_2N_array;			// pair-nucleon state angular momentum
	int* S_2N_array;			// pair-nucleon state total spin
	int* J_2N_array;			// pair-nucleon state total angular momentum
	int* T_2N_array;			// pair-nucleon state total isospin
	int* L_1N_array;			// orbital-nucleon state angular momentum
	int* two_J_1N_array;		// orbital-nucleon state total angular momentum x2
	int* two_J_3N_array;		// three-nucleon state total angular momentum x2
	int* two_T_3N_array;		// three-nucleon state total isospin x2
	int* P_3N_array;			// three-nucleon state parity
	int* chn_3N_idx_array;		// Indices of each 3N JP-channel
	int  N_chn_3N;				// Number of distinct 3N JP-channels
} pw_3N_statespace;

/* Free wave-packet (FWP) statespace struct */
typedef struct fwp_statespace{
	int 	Np_WP;				// Number of p-momentum WPs
	int 	Nq_WP;				// Number of q-momentum WPs
	int 	Np_per_WP;			// Number of p-momentum quadrature points in each p-momentum WP
	int 	Nq_per_WP;			// Number of q-momentum quadrature points in each q-momentum WP
	double* p_WP_array;			// Boundaries of p-momentum WPs
	double* q_WP_array;			// Boundaries of q-momentum WPs
	double* p_array;			// p-momentum quadrature points/average momenta, for all bins
	double* wp_array;			// p-momentum quadrature weights, for all bins
	double* fp_array;			// p-momentum weight-function at p_array-points, for all bins
	double* q_array;			// q-momentum quadrature points/average momenta, for all bins
	double* wq_array;			// q-momentum quadrature weights, for all bins
	double* fq_array;			// q-momentum weight-function at q_array-points, for all bins
	double* norm_p_array;		// p-momentum normaliszation, for all bins
	double* norm_q_array;		// q-momentum normaliszation, for all bins
} fwp_statespace;

/* Scattering wave-packet (SWP) statespace struct */
typedef struct swp_statespace{
	int 	Np_WP;				// Number of p-momentum WPs
	int 	Nq_WP;				// Number of q-momentum WPs
	int 	num_2N_unco_states;	// Number of uncoupled NN states
	int 	num_2N_coup_states;	// Number of coupled NN states
	double  E_bound;			// Deuteron bound state energy (with negative sign)
	double* e_SWP_unco_array;	// Boundaries of uncoupled scattering WPs
	double* e_SWP_coup_array;	// Boundaries of coupled scattering WPs
	double* C_SWP_unco_array;	// Basis transformation matrices for uncoupled WPs
	double* C_SWP_coup_array;	// Basis transformation matrices for coupled WPs
	double* q_WP_array;			// Boundaries of q-momentum WPs
} swp_statespace;

/* Struct containing all information regarding on-shell indices for on-shell
 * energies and corresponding bins, as well as deuteron channels, for ALL 3N channels */ 
typedef struct solution_configuration{
	size_t  num_T_lab;				// Number of on-shell bins/energies to calculate
	double* T_lab_array;			// On-shell lab  energies  (T_lab)
	double* q_com_array;			// On-shell c.m. q-momenta (q_com)
	double* E_com_array;			// On-shell c.m. energies  (E_com)
	int*    q_com_idx_array;		// Index-array for on-shell q-bin (used for elastic scattering)
	int*    alphapq_idx_array;		// Index-array for on-shell alpha-,p-,q-bins (used for breakup)
	int**   deuteron_idx_arrays;	// Index-arrays deuteron-channels in all 3N-channels
	int*    deuteron_num_array;		// Contains number of deuteron-channels in all 3N-channels
} solution_configuration;

/* Struct containing all information regarding on-shell indices for on-shell
 * bins, as well as deuteron channels, for GIVEN 3N channel */ 
typedef struct channel_os_indexing{
	size_t  num_T_lab;				// Number of elastic on-shell bins/energies to calculate
	int*    q_com_idx_array;		// Index-array for elastic on-shell q-bins
	int*    deuteron_idx_array;		// Index-array for deuteron-channels in given 3N-channel
	int     num_deuteron_states;	// Number of deuteron-channels in given 3N-channel
	size_t* OS_BU_idx_array;		// on-shell breakup index array
	int*    pq_com_idx_array;		// Index-array for all on-shell p- and q-bins (elastic AND breakup)
} channel_os_indexing;

typedef struct run_params{
	int         two_J_3N_max;
	int         Np_WP;
	int         Nq_WP;
	int         J_2N_max;
	int         Nphi;
	int         Nx;
	int         Np_per_WP;
	int         Nq_per_WP;
	int			channel_idx;
	int			P123_omp_num_threads;
	int			max_TFC;
	double 		chebyshev_t;
	double 		chebyshev_s;
	bool        parallel_run;
	bool		P123_recovery;
	bool 		tensor_force;
	bool		isospin_breaking_1S0;
	bool 		midpoint_approx;
	bool		calculate_and_store_P123;
	bool		solve_faddeev;
	bool		production_run;
	bool	    parameter_walk;
	int			PSI_start;
	int			PSI_end;
	std::string potential_model;
	std::string subfolder;
	std::string p_grid_type;
	std::string p_grid_filename;
	std::string q_grid_type;
	std::string q_grid_filename;
	std::string parameter_file;
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

