program convert

implicit none

integer :: i,j,k
integer :: nsteps, nspec, maxatoms
integer :: idu
integer, allocatable, dimension(:) :: natoms

double precision, allocatable, dimension(:,:) :: x,y,z

character*2, allocatable, dimension(:,:) :: symb
character*80 :: lammpsfile,positionfile

open(10,file='convert.inpt')
read(10,*)nsteps
read(10,*)nspec
allocate(natoms(nspec))
maxatoms=0
do i=1,nspec
   read(10,*)natoms(i)
   if(natoms(i).gt.maxatoms)maxatoms=natoms(i)
enddo
allocate(x(nspec,maxatoms),y(nspec,maxatoms),z(nspec,maxatoms))
allocate(symb(nspec,maxatoms))
read(10,'(a)')lammpsfile
open(20,file=lammpsfile)
read(10,'(a)')positionfile
open(30,file=positionfile)

do i=1,nsteps

   do j=1,9
      read(20,*)
   enddo   
   do j=1,nspec
      do k=1,natoms(j)
         read(20,*)symb(j,k),idu,idu,x(j,k),y(j,k),z(j,k)
      enddo
   enddo   

   do j=1,nspec
      do k=1,natoms(j)
         if(symb(j,k).eq.'Cl')write(30,*)x(j,k),y(j,k),z(j,k)
      enddo
   enddo   
   do j=1,nspec
      do k=1,natoms(j)
         if(symb(j,k).eq.'La')write(30,*)x(j,k),y(j,k),z(j,k)
      enddo
   enddo   
   do j=1,nspec
      do k=1,natoms(j)
         if(symb(j,k).eq.'Mg')write(30,*)x(j,k),y(j,k),z(j,k)
      enddo
   enddo   
   do j=1,nspec
      do k=1,natoms(j)
         if(symb(j,k).eq.'Na')write(30,*)x(j,k),y(j,k),z(j,k)
      enddo
   enddo   
enddo
close(20)
close(30)

end program

