Installation Instructions:
	
- Prerequisites:
	> Fortran compiler (e.g., gfortran)

- Steps to compile and run:
 	> Ensure working at correct directory 
	    mv /src
	> Compile the program:
            gfortran md.f90 -o md.x
	> Run the program:
	    ./md.x
	> Check energy and trajectory files:
 	    vim energy.out
	    vim traj.xyz
	> Plot energy with your desired software, e.g.:
	    gnuplot energy.out
	> Visualize trajectory with your desired software, e.g:
  	    vmd traj.xyz
	> Make sure to change parameters (Nsteps,sigma,epsilon) and have fun visualizing particle motion =)

