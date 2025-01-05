Project 1: Hartree-Fock and Møller-Plesset pertubation theory order 2 (MP2) methods implementation
==================================================================================================
*Nicole Potenza, Davi Augusto Matosso, Miguel Angel Díaz Gutiérrez*

The following code runs a Hartree-Fock calculations followed by an MP2 correction for the energy. This is achieved by using the following equations:

$$
\begin{equation}
E_{HF}=E_{NN}+2\sum_{i=1}^{N_{occ}}\langle i | h | i\rangle + \sum_{i=1}^{N_occ}\sum_{j=1}^{N_{occ}} [2\langle ij|ij\rangle-\langle ij|ji\rangle]
\end{equation}
$$

where $E_{NN}$ is the nuclear repulsion, the first sum refers to the one electron integrals while the terms in the squared brackets refer to the tow-electron integrals
