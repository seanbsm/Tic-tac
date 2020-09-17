
#include "read_psi.h"

using namespace std;

void read_psi (char* filename, double** psi_read, int *Np, double **p, double **wp, int *Nq, double **q, double **wq, int *Nalpha, int **get_L, int **get_S, int **get_J, int **get_T, int **get_l, int **get_two_j){

    bool print_content = false;

    FILE *pFile;
    cout << "Read file: " << filename << "\n";
    
    pFile = fopen (filename, "r");
    if (fscanf(pFile, "%d %d %d", Np, Nq, Nalpha) != 3){ raise_error("problem_reading_data"); }

    if (print_content) cout << "Np = " << *Np << ", Nq = " << *Nq << ", Nalpha = " << *Nalpha << "\n";

    double plocal[*Np];
    double wplocal[*Np];
    double qlocal[*Nq];
    double wqlocal[*Nq];
    double psilocal[(*Np) * (*Nq) * (*Nalpha)];

    int get_L_local[*Nalpha];
    int get_S_local[*Nalpha];
    int get_J_local[*Nalpha];
    int get_T_local[*Nalpha];
    int get_l_local[*Nalpha];
    int get_two_j_local[*Nalpha];

    (*p) = new double[*Np];
    (*wp) = new double[*Np];

    (*q) = new double[*Nq];
    (*wq) = new double[*Nq];

    *get_L = new int[*Nalpha];
    *get_S = new int[*Nalpha];
    *get_J = new int[*Nalpha];
    *get_T = new int[*Nalpha];
    *get_l = new int[*Nalpha];
    *get_two_j = new int[*Nalpha];

    *psi_read = new double[(*Np) * (*Nq) * (*Nalpha)];

    for (int p_index = 0; p_index <= *Np - 1; p_index++){
        if (fscanf(pFile, "%lf %lf", &plocal[p_index], &wplocal[p_index]) != 2){ raise_error("problem_reading_data"); }
        (*p)[p_index] = plocal[p_index];
        (*wp)[p_index] = wplocal[p_index];

        if (print_content) cout << plocal[p_index] << " " << wplocal[p_index] << " " << *Np << "\n";
    }

    for (int q_index = 0; q_index <= *Nq - 1; q_index++){
        if (fscanf(pFile, "%lf %lf", &qlocal[q_index], &wqlocal[q_index]) != 2){ raise_error("problem_reading_data"); }
        (*q)[q_index] = qlocal[q_index];
        (*wq)[q_index] = wqlocal[q_index];

        if (print_content) cout << qlocal[q_index] << " " << wqlocal[q_index] << "\n";
    }

    if (print_content) cout << "L12 S12 J12 T12 l3 2j3\n";
    
    for (int alpha = 0; alpha <= *Nalpha - 1; alpha++){
        if (fscanf(pFile, "%d %d %d %d %d %d\n", &get_L_local[alpha], &get_S_local[alpha], &get_J_local[alpha], &get_T_local[alpha], &get_l_local[alpha], &get_two_j_local[alpha]) != 6){ raise_error("problem_reading_data"); }
        (*get_L)[alpha] = get_L_local[alpha];
        (*get_S)[alpha] = get_S_local[alpha];
        (*get_J)[alpha] = get_J_local[alpha]; 
        (*get_T)[alpha] = get_T_local[alpha];
        (*get_l)[alpha] = get_l_local[alpha];
        (*get_two_j)[alpha] = get_two_j_local[alpha];

        if (print_content) cout << get_L_local[alpha] << "   " << get_S_local[alpha] << "   " << get_J_local[alpha] << "   " << get_T_local[alpha] << "   " << get_l_local[alpha] << "   " << get_two_j_local[alpha] << "\n";
    }

    for (int p_index = 0; p_index <= (*Np) - 1; p_index++){
        for (int q_index = 0; q_index <= (*Nq) - 1; q_index++){
            for (int alpha = 0; alpha <= (*Nalpha) - 1; alpha++){
                if (fscanf(pFile, "%lf\n", &psilocal[p_index * (*Nq) * (*Nalpha) + q_index * (*Nalpha) + alpha]) != 1){ raise_error("problem_reading_data"); }
                (*psi_read)[p_index * (*Nq) * (*Nalpha) + q_index * (*Nalpha) + alpha] = psilocal[p_index * (*Nq) * (*Nalpha) + q_index * (*Nalpha) + alpha];
            }
        }
    }
    fclose(pFile);
}

void get_all_states(double** state_3N_symm_array, double** state_3N_asym_array, int& Np, double** p, double** wp, int& Nq, double** q, double** wq, int& Nalpha, int** L_2N, int** S_2N, int** J_2N, int** T_2N, int** l_3N, int** two_j_3N){

    /* Read wave functions to argument-arrays */
    char psi_3H_filename[100];
    sprintf(psi_3H_filename, "%s", "Triton_states/H3_psi_N3LO_EM500.dat");
    read_psi(psi_3H_filename, state_3N_symm_array, &Np, p, wp, &Nq, q, wq, &Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N);
    
    char psiasymm_3H_filename[100];
    sprintf(psiasymm_3H_filename, "%s", "Triton_states/H3_psiasymm_N3LO_EM500.dat");
    read_psi(psiasymm_3H_filename, state_3N_asym_array, &Np, p, wp, &Nq, q, wq, &Nalpha, L_2N, S_2N, J_2N, T_2N, l_3N, two_j_3N);
}
