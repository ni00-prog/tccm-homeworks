!MOLECULAR DYNAMICS SIMULATION OF ARGON ATOMS USING LENNARD-JONES POTENTIAL
program md
implicit none
real*8,allocatable :: coord(:,:),mass(:),dist(:,:),vel(:,:),acc(:,:),dtcoord(:,:),old_acc(:,:)
real*8 :: r,V,Vtot,K,Ktot,E,Etot
real*8,parameter :: epsilon=0.0661,sigma=0.3345,dt=0.2
integer :: n,Natoms,read_Natoms,i,j 
integer,parameter :: Nsteps=1000,Ntraj=50,Nenergy=5  ! Number of steps of MD ; Output traj/energy every x steps;
character(len=30) :: input_file,output_file,energy_file

! Get number of atoms
input_file='inp.txt'
Natoms=read_Natoms(input_file)

! Read file (coordinates and masses) 
allocate(coord(Natoms,3),mass(Natoms))
call read_molecule(input_file,Natoms,coord,mass)

! Calculate distances
allocate(dist(Natoms,Natoms))
call compute_distances(Natoms,coord,dist)  

! Print initial conditions (Natoms,coord,masses,distances)
write(*,*) '##############################'
write(*,*) '#### INITIAL CONDITIONS ######'
write(*,*) '##############################'
write(*,*) 'Natoms: ', Natoms
write(*,*) 'Nsteps: ', Nsteps
write (*,*) 'Coordinates:'
do i=1,Natoms
      write(*,130) (coord(i,j), j=1,3)
enddo
write(*,*) 'Masses: '
write(*,130) mass
write (*,*) 'Distances:'
do i=1,Natoms
   write(*,130) (dist(i,j), j=1,Natoms)
enddo

! Calculate potential 

Vtot= V(epsilon,sigma,dist,Natoms)
write(*,110) 'Potential Energy: ',Vtot

! Calculate kinetic energy

allocate(vel(Natoms,3))
vel=0.d0
Ktot= K(Natoms,vel,mass)
write(*,110) 'Kinetic Energy: ', Ktot

! Calculate total energy

Etot= Ktot + Vtot   !direct calculation instead of function call
write(*,110) 'Total Energy: ', Etot

! Calculate acceleration

allocate(acc(Natoms,3))
call compute_acc(Natoms,coord,mass,dist,acc)
write(*,*) 'Acceleration: '

do i=1,Natoms
    write(*,130) (acc(i,j), j=1,3)
enddo

! Create output files

output_file='traj.xyz'
open(13, file=output_file,status='replace') !create trajectories's output file
close(13)

energy_file='energies.out'
open(14,file=energy_file, action='write', status='replace')     !create energies's output file
write(14,*) 'Step    kineticE   potentialE   TotalE'
close(14)

write(*,*) 'Trajectory are found in: ',output_file
write(*,*) 'Energies are found in: ',energy_file

!Write initial configuration
open(13, file=output_file,position='append')
write(13,*) Natoms
write(13,120) 'Step 0, Total E=',Etot
do i=1,Natoms
    write(13,140) 'Ar',(coord(i,j), j=1,3)
enddo
close(13)

!Main MD propagation loop 

allocate(old_acc(Natoms,3))
old_acc = acc

do n=1,Nsteps
    old_acc = acc !store old acceleration

!Update position Verlet algorithm eqz.6
    dtcoord = vel*dt + 0.5d0*acc*(dt**2)
    coord = coord + dtcoord                       !compute new acceleration
    
    call compute_distances(Natoms,coord,dist)     !recompute distances with new positiond
    call compute_acc(Natoms,coord,mass,dist,acc)  !compute new acceleration

    !Update velocities with Verlet eqz.7
    vel=vel+0.5d0*(old_acc + acc)*dt

    !calculate and save energies every Nenergy steps
    if (mod(n,Nenergy) == 0) then
       Ktot = K(Natoms,vel,mass)
       Vtot = V(epsilon,sigma,dist,Natoms)
       Etot = Ktot + Vtot
       call write_energies(energy_file,n,Ktot,Vtot,Etot)
    endif

    !save trajectory every Ntraj steps
    if (mod(n,Ntraj) == 0) then
       Ktot = K(Natoms,vel,mass)                          
       Vtot = V(epsilon,sigma,dist,Natoms)                
       Etot = E(Ktot,Vtot)                      
       call write_traj(n,Natoms,coord,output_file,Etot)
    endif
enddo

deallocate(coord,mass,dist,vel,acc,old_acc,dtcoord)
100 format (A,3f12.6)
110 format(A,F12.6)
120 format (A,F15.9)
130 format(3(2x,f12.6))    
140 format (A,3(f15.9))
end program md

! Obtain Natoms 
integer function read_Natoms(input_file) result(Natoms)
character(len=20),intent(in) :: input_file
open(10,file=input_file)
read(10,*) Natoms
close(10)
end function read_Natoms

! Obtain coordinates and mass
subroutine read_molecule(input_file, Natoms, coord, mass)
implicit none
character(len=20), intent(in) :: input_file
integer, intent(in) :: Natoms
real*8, intent(out) :: coord(Natoms,3)
real*8, intent(out) :: mass(Natoms)
integer :: i,j

! Read data
open(11,file=input_file)
read(11,*)
do i=1,Natoms
   read(11,*) coord(i,:),mass(i)
enddo
close(11)

end subroutine read_molecule

! Compute distances
subroutine compute_distances(Natoms, coord, dist)
implicit none
integer, intent(in) :: Natoms
real*8, intent(in) :: coord(Natoms,3)
real*8, intent(out) :: dist(Natoms,Natoms)
integer :: i,j,k
do i=1,Natoms
   do j=1,Natoms
      dist(i, j) = sqrt((coord(i, 1) - coord(j, 1))**2 +(coord(i, 2) - coord(j, 2))**2 +(coord(i, 3) - coord(j, 3))**2)
   enddo
enddo

end subroutine compute_distances
! Calculate potential energy
real*8 function V(epsilon, sigma, dist, Natoms)
    implicit none
    integer, intent(in) :: Natoms
    real*8, intent(in) :: epsilon, sigma, dist(Natoms,Natoms)
    real*8 :: r
    integer :: i, j
    
    V = 0.d0
    do i = 1, Natoms
        do j = i+1, Natoms
            r = dist(i,j)
            if (r <= 0.d0) then
                V = 1.0d10
            else
                V = V + 4*epsilon*((sigma/r)**12 - (sigma/r)**6)  !Eq. 2
            endif
        enddo
    enddo
    return
end function V

! Calculate kinetic energy
real*8 function K(Natoms, vel, mass)
    implicit none
    integer, intent(in) :: Natoms
    real*8, intent(in) :: vel(Natoms,3)
    real*8, intent(in) :: mass(Natoms)
    integer :: i, j
    
    K = 0.d0
    do i = 1, Natoms
        do j = 1, 3
            K = K + 0.5d0*mass(i)*(vel(i,j)**2)  ! Eq. 3
        enddo
    enddo
    return
end function K

! Calculate total energy
real*8 function E(Ktot, Vtot)
    implicit none
    real*8, intent(in) :: Ktot, Vtot
    E = Ktot + Vtot  !Eq 4
    return
end function E

! Calculate acceleration
subroutine compute_acc(Natoms, coord, mass, dist, acc)
    implicit none
    integer, intent(in) :: Natoms
    real*8, intent(in) :: coord(Natoms,3)
    real*8, intent(in) :: mass(Natoms)
    real*8, intent(in) :: dist(Natoms,Natoms)
    real*8, intent(out) :: acc(Natoms,3)
    real*8, parameter :: epsilon = 0.0661, sigma = 0.3345
    real*8 :: r(3), rij, U
    integer :: i, j, k
    
    acc = 0.d0
    do i = 1, Natoms
        do j = i+1, Natoms
            r = coord(i,:) - coord(j,:)
            rij = dist(i,j)
            if (rij > 1.d-10) then
                U = 24.d0 * epsilon / rij * ((sigma/rij)**6 - 2*(sigma/rij)**12)  !Eq. 6
                do k = 1, 3
                    acc(i,k) = acc(i,k) - (1.d0/mass(i))*U*r(k)/rij  !Eq 5 atom i
                    acc(j,k) = acc(j,k) + (1.d0/mass(j))*U*r(k)/rij  !Opposite force for atom j
                enddo
            endif
        enddo
    enddo
end subroutine compute_acc

subroutine write_energies(energy_file,step,Ktot,Vtot,Etot)
implicit none
integer,intent(in)              :: step
real*8,intent(in)               :: Ktot,Vtot,Etot
character(len=30),intent(inout) :: energy_file

open(14,file=energy_file,action='write',position='append')
write(14,100) step,Ktot,Vtot,Etot
close(14)

100 format (I5,3f15.9)
end subroutine write_energies

subroutine write_traj(step,Natoms,coord,output_file,Etot)
implicit none
integer,intent(in)              :: step,Natoms
real*8,intent(in)               :: coord(Natoms,3),Etot
character(len=30),intent(inout) :: output_file
integer                         :: i,j

open(13,file=output_file,action='write',position='append')
write(13,100) Natoms
write(13, 110) 'Step:', step, '; Total E =', Etot
do i = 1, Natoms
   write(13, 120) 'Ar', (coord(i,j),j=1,3)
end do

close(13)
100 format (I3)
110 format(A,I5,A,f15.9)
120 format (A,3(f15.9))
end subroutine write_traj

