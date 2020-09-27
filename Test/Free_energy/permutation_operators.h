#ifndef PERMUTATION_OPERATORS_H
#define PERMUTATION_OPERATORS_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <vector>

/* Modules used in benchmark code */
#include"mkl.h"
#include<math.h>
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

typedef struct pw_table
{
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

typedef struct alpha_table
{
    int alpha_alphaprime_index;
    int alpha_index;
    int alphaprime_index;
} alpha_table;


typedef struct pmesh_table
{
    int index_ptable;
    double mesh_ptable;
    double weight_ptable;
} pmesh_table;


typedef struct qmesh_table
{
    int index_qtable;
    double mesh_qtable;
    double weight_qtable;
} qmesh_table;

void read_P123_h5_data_file(double* P123);
void read_P123_csv_data_file(double* P123);
void read_P123_bin_data_file(double* P123, int P123_array_size);

void calculate_antisymmetrization_operator(int &Np, int &Nq, int& Nalpha,
                                           double** A123);

#endif // PERMUTATION_OPERATORS_H