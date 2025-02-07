// project 1 - hf and MP2 energy calculation 
#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>
#include <stdint.h>

// function to find an integral in the sparse storage
double find_integral(int i, int j, int k, int l, int64_t n_integrals, 
                    int32_t* index, double* value) {
   	for (int64_t n = 0; n < n_integrals; n++) {		//check all permutations
        int idx_i = index[4*n];
        int idx_j = index[4*n+1];
        int idx_k = index[4*n+2];
        int idx_l = index[4*n+3];
        
        if ((idx_i == i && idx_j == j && idx_k == k && idx_l == l) ||  //check equivalent permutations
            (idx_i == i && idx_j == l && idx_k == k && idx_l == j) ||
            (idx_i == k && idx_j == j && idx_k == i && idx_l == l) ||
            (idx_i == k && idx_j == l && idx_k == i && idx_l == j) ||
            (idx_i == j && idx_j == i && idx_k == l && idx_l == k) ||
            (idx_i == j && idx_j == k && idx_k == l && idx_l == i) ||
            (idx_i == l && idx_j == i && idx_k == j && idx_l == k) ||
            (idx_i == l && idx_j == k && idx_k == j && idx_l == i)) {
            return value[n];
        }
    }
    return 0.0; 
}
//Hartree-Fock energy calculation
double calculate_hf_energy(trexio_t* file) {
    trexio_exit_code rc;
    double nuclear_repulsion;
    int32_t mo_num, n_up;
    double *one_e_ints;
    
    //read nuclear repulsion
    rc = trexio_read_nucleus_repulsion(file, &nuclear_repulsion);
    printf("Nuclear repulsion energy: %.7f\n", nuclear_repulsion);
    
    //number of electrons 
    rc = trexio_read_electron_up_num(file, &n_up);
    int n_electrons = 2 * n_up; 				 // Total n. electrons for closed shell
    
    //number of MOs
    rc = trexio_read_mo_num(file, &mo_num);
    printf("Number of electrons: %d\n", n_electrons);
    printf("Number of MOs: %d\n", mo_num);
    
    //read one-electron integrals
    one_e_ints = malloc(mo_num * mo_num * sizeof(double));
    if (one_e_ints == NULL) {
        fprintf(stderr, "Malloc failed for one-electron integrals\n");
        exit(1);
    }
    rc = trexio_read_mo_1e_int_core_hamiltonian(file, one_e_ints);

    //read two-electron integrals
    int64_t n_integrals;
    rc = trexio_read_mo_2e_int_eri_size(file, &n_integrals);
    printf("Number of two-electron integrals: %lld\n", n_integrals);
    
    int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
    double* value = malloc(n_integrals * sizeof(double));
    if (index == NULL || value == NULL) {
        fprintf(stderr, "Malloc failed for two-electron integrals\n");
        free(one_e_ints);
        free(index);
        free(value);
        exit(1);
    }
    
    int64_t buffer_size = n_integrals;
    rc = trexio_read_mo_2e_int_eri(file, 0, &buffer_size, index, value);
    
    //calculate energy components
    double energy = nuclear_repulsion;
    double one_e_contribution = 0.0;
    
    //one-electron contribution 
    for (int i = 0; i < n_up; i++) {
        one_e_contribution += 2.0 * one_e_ints[i * mo_num + i];   //factor 2 for closed shell
    }
    energy += one_e_contribution;
    printf("One-electron contribution: %.7f\n", one_e_contribution);
    
    //two-electron contribution
    double two_e_contribution = 0.0;
    for (int i = 0; i < n_up; i++) {
        for (int j = 0; j < n_up; j++) {
            double Jij = find_integral(i, j, i, j, n_integrals, index, value);
            double Kij = find_integral(i, j, j, i, n_integrals, index, value);
            two_e_contribution += 2.0 * Jij - Kij;
        }
    }
    energy += two_e_contribution;
    printf("Two-electron contribution: %.7f\n", two_e_contribution);
    
    free(one_e_ints);
    free(index);
    free(value);
    
    return energy;
}  
//second order Moller-Plesset perturbation correction
double calculate_mp2_correction(trexio_t* file){
  trexio_exit_code rc;
  int32_t mo_num, n_up;
  // read the number of electrons and number of orbitals
  rc = trexio_read_mo_num(file, &mo_num);
  rc = trexio_read_electron_up_num(file, &n_up);
  // read the orbital energies
  double *eps = malloc(mo_num * sizeof(double));
  rc = trexio_read_mo_energy(file, eps);
  if (rc != TREXIO_SUCCESS) {
      fprintf(stderr, "The energies could not be read\n");
      free(eps);
      return 0.0;
    } //if reading fails.
// read the two-electron integrals
  int64_t n_integrals;
  rc = trexio_read_mo_2e_int_eri_size(file, &n_integrals);

  int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
  double* value = malloc(n_integrals * sizeof(double));
  int64_t buffer_size = n_integrals;
  rc = trexio_read_mo_2e_int_eri(file, 0, &buffer_size, index, value);

// computing of the mp2 energy (using nested loops)
  double mp2_energy = 0.0;
  // first bucle over the occupied orbitals
  for(int i = 0; i < n_up; i++) {
    for(int j = 0; j < n_up; j++) {
      // second bucle over the virtual orbitals (k, l)
      for(int a = n_up; a < mo_num; a++) {
        for (int b = n_up; b < mo_num; b++) {
            // orbital energies
            double denom = eps[i] + eps[j] - eps[a] - eps[b];
            // two electrons integrals
            double ijab = find_integral(i, j, a, b, n_integrals, index, value);
            double ijba = find_integral(i, j, b, a, n_integrals, index, value);
            // the energy following the equation
            mp2_energy += (ijab * (2.0*ijab -ijba)) / denom;
          }
        }
      }
    }
    free(eps);
    free(index);
    free(value);

  return mp2_energy;
}
//calculate both HF and MP2 energies for H2O molecule
int main() {
  // open the file
    trexio_exit_code rc;
    trexio_t* file = trexio_open("../data/h2o.h5", 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Could not open h2o.h5\n");
        return 1;
    }
    // Calling the functions and printing the results. 
    double hf_energy = calculate_hf_energy(file);
    printf("\nH2O Hartree-Fock Total Energy: %.7f atomic units\n", hf_energy);

    double mp2_correction = calculate_mp2_correction(file);
    printf("MP2 Energy correction: %.7f atomic units\n", mp2_correction);
    printf("Total energy (HF + MP2): %.7f atomic units\n", hf_energy + mp2_correction);
  
    trexio_close(file);
    return 0;
}
