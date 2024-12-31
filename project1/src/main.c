#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>
#include <stdint.h>

double calculate_hf_energy(trexio_t* file) {
    trexio_exit_code rc;
    double nuclear_repulsion;
    int32_t mo_num, n_electrons;
    double *one_e_ints;
    
    rc = trexio_read_nucleus_repulsion(file, &nuclear_repulsion);
    printf("Nuclear repulsion energy: %.7f\n", nuclear_repulsion);
    
    rc = trexio_read_mo_num(file, &mo_num);
    rc = trexio_read_electron_up_num(file, &n_electrons);
    printf("Number of electrons: %d\n", n_electrons);
    printf("Number of MOs: %d\n", mo_num);
    
    one_e_ints = malloc(mo_num * mo_num * sizeof(double));
    if (one_e_ints == NULL) {
        fprintf(stderr, "Malloc failed for one-electron integrals\n");
        exit(1);
    }
    rc = trexio_read_mo_1e_int_core_hamiltonian(file, one_e_ints);

    int64_t n_integrals;
    rc = trexio_read_mo_2e_int_eri_size(file, &n_integrals);
    printf("Number of two-electron integrals: %lld\n", n_integrals);
    
    int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
    if (index == NULL) {
        fprintf(stderr, "Malloc failed for index\n");
        free(one_e_ints);
        exit(1);
    }
    
    double* value = malloc(n_integrals * sizeof(double));
    if (value == NULL) {
        fprintf(stderr, "Malloc failed for value\n");
        free(one_e_ints);
        free(index);
        exit(1);
    }
    
    int64_t buffer_size = n_integrals;
    rc = trexio_read_mo_2e_int_eri(file, 0, &buffer_size, index, value);
    
    double energy = nuclear_repulsion;
    double one_e_contribution = 0.0;
    
    for (int i = 0; i < n_electrons; i++) {
        one_e_contribution += 2.0 * one_e_ints[i * mo_num + i];
    }
    energy += one_e_contribution;
    printf("One-electron contribution: %.7f\n", one_e_contribution);
    
    double two_e_contribution = 0.0;
    for (int i = 0; i < n_electrons; i++) {
        for (int j = 0; j < n_electrons; j++) {
            for (int64_t n = 0; n < buffer_size; n++) {
                int k = index[4*n+2];
                int l = index[4*n+3];
                
                if (index[4*n] == i && index[4*n+1] == j) {
                    if (k == i && l == j) {
                        two_e_contribution += 2.0 * value[n];
                    }
                    if (k == j && l == i) {
                        two_e_contribution -= value[n];
                    }
                }
            }
        }
    }
    energy += two_e_contribution;
    printf("Two-electron contribution: %.7f\n", two_e_contribution);
    
    free(one_e_ints);
    free(index);
    free(value);
    
    return energy;
}

int main() {
    trexio_exit_code rc;
    trexio_t* file = trexio_open("../data/h2o.h5", 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Could not open h2o.h5\n");
        return 1;
    }
    
    double hf_energy = calculate_hf_energy(file);
    printf("\nH2O Hartree-Fock Total Energy: %.7f atomic units\n", hf_energy);
    
    trexio_close(file);
    return 0;
}

//the hf energy is wrong. 
