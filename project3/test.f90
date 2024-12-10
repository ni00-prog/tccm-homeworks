PROGRAM main_md									!name program without space
    implicit none

    ! Variable declarations
    real*8, allocatable :: coord(:,:), mass(:), dist(:,:), vel(:,:), acc(:,:), dtcoord(:,:)
    real*8 :: r, V, Vtot, K, Ktot, E, Etot
    real*8, parameter :: epsilon=0.0661, sigma=0.3345, dt=0.2
    integer :: n, Natoms, read_Natoms, i, j
    integer, parameter :: Nsteps=10, Ntraj=10, Nenergy=1 			!number of steps of MD; output trj/energy every x steps
    character(len=30) :: input_file, output_file, energy_file

    !INITIALIZATION
    
    ! Get number of atoms from input file
    input_file= 'inp.txt'
    Natoms=read_Natoms(input_file)						!read_Natoms is a function not a variable

    !Read file (coordinate and masses)
    		  
    allocate(dist(Natoms,Natoms))						!add other allocates

    call compute_distances(Natoms,coord,dist)					!add other calls 

    !PRINT INITIAL CONDITIONS

    write(*,*) 'INITIAL CONDITIONS'
    write(*,*) 'Natoms: ',Natoms
    write(*,*) 'Nsteps: ',Nsteps
    write(*,*) 'Coordinates (x,y,z):'
    
       do i=1, Natoms
          write(*,130) (coord(i,j), j=1,Natoms)
       end do

    write(*,*) 'Masses: '
    write(*,130) mass
    write(*,*) 'Distances:'
       
       do i=1, Natoms
          write(*,130) (dist(i,j), j=1,Natoms)
       enddo
    
130 format(3F12.6) 								!add the format

END PROGRAM main_md
