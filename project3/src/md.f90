program md
implicit none
real*8,allocatable :: coord(:,:),mass(:),dist(:,:),vel(:,:),acc(:,:),dtcoord(:,:)
real*8 :: r,V,Vtot,K,Ktot,E,Etot
real*8,parameter :: epsilon=0.0661,sigma=0.3345,dt=0.2
integer :: n,Natoms,read_Natoms,i,j 
integer,parameter :: Nsteps=10,Ntraj=10,Nenergy=1  ! Number of steps of MD ; Output traj/energy every x steps;
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
      write(*,130) (coord(i,j), j=1,Natoms)
enddo
write(*,*) 'Masses: '
write(*,130) mass
write (*,*) 'Distances:'
do i=1,Natoms
   write(*,130) (dist(i,j), j=1,Natoms)
enddo

130 format(3(2x,f12.6))

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



