//libraries
#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>

//HF energy
double calculate_hf_energy(trexio_t* file) {
    double nuclear_repulsion;			//nuclear repulsion energy
    int32_t mo_num;				//number MOs
    int32_t n_electrons;			//number electrons
    double *one_e_ints;				//one-elec. integrals array

//read nuclear repulsion energy from h2o.h5
    rc = trexio_read_nucleus_repulsion(file, &nuclear_repulsion);
    
//read MOs and electrons
    rc = trexio_read_mo_num(file, &mo_num);
    rc = trexio_read_electron_up_num(file, &n_electrons);
    
//read one electron integrals
   one_e_ints = malloc(mo_num * mo_num * sizeof(double));		//allocate memory
   rc = trexio_read_mo_1e_int_core_hamiltonian(file, one_e_ints);

//read two electrons integrals
   


