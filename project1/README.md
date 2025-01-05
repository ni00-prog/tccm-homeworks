Project 1: Hartree-Fock and Møller-Plesset pertubation theory order 2 (MP2) methods implementation
==================================================================================================
*Nicole Potenza, Davi Augusto Matosso, Miguel Angel Díaz Gutiérrez*

The following code runs a Hartree-Fock calculations followed by an MP2 correction for the energy. This is achieved by using the following equations:

$$
\begin{equation}
E_{HF}=E_{NN}+2\sum_{i=1}^{N_{occ}}\langle i | h | i\rangle + \sum_{i=1}^{N_occ}\sum_{j=1}^{N_{occ}} [2\langle ij|ij\rangle-\langle ij|ji\rangle]
\end{equation}
$$

where $E_{NN}$ is the nuclear repulsion, the first sum refers to the one electron integrals while the terms in the squared brackets refer to the two-electron integrals. The MP2 correction has the form:

$$
\begin{equation}
E_{MP2}=\sum_{(i,j)\in occupied}\sum_{a,b\in virtuals}\dfrac{\langle ij|ab\rangle(2\langle ij|ab \rangle - \langle ij|ba\rangle)}{e_i + e_j - e_a -e_b}
\end{equation}
$$

where $e_i$, $e_j$ are the occupied orbital energies and $e_a$, $e_b$ are the virtual orbital energies. 

This code uses TREXIO to extract quantum chemical information such as: 

* One-electron and two-electrons integrals.
* Number of electrons.
* Number of molecular orbitals.
* Number of occupied orbitals.
* Nuclear repulsion energies.
* Orbital energies.

Which are used for both methods. 

Structure of the directory
---------------------------
The main directory "Project 1" is part of the repository tccm-homeworks. 
