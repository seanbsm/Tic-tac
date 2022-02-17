#ifndef ERROR_MANAGEMENT_H
#define ERROR_MANAGEMENT_H

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

/* hdf5-storage libaries */
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

/* MKL libraries */
#include "mkl.h"
#include "mkl_spblas.h"

void raise_error (std::string error_cause);
void check_h5_read_call(herr_t ret);
void check_h5_read_table_call(herr_t ret);
void check_h5_close_call(herr_t ret);
bool verify_sparse_mkl_success(sparse_status_t mkl_sparse_status, bool print_true);

#endif // ERROR_MANAGEMENT_H