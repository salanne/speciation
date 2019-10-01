program polym

! analyze polymer species in mixtures

implicit none

integer :: nmax
parameter (nmax=800)

integer :: i,j,k,l,ll
integer :: imark,jmark,nskip
integer :: nconfigs,nequil,nanion1,ncationpol1,ncationfree1,nion,nanion2
integer :: ncationpol2
integer :: nneighb,nshar,nchain,nat,ichain,jchain,ntri,nmol
integer, allocatable, dimension(:,:,:,:) :: spec
integer, allocatable, dimension(:) :: nneighbcat2

integer, dimension(10)        :: shar
integer, dimension(nmax)      :: chain,nan1,nan2,ncat1,ncat2,ncharge
integer, dimension(nmax,nmax) :: atoms

double precision :: boxlength,halfbox,halfboxrec
double precision :: du,dx,dy,dz,r2,totfrac1,totfrac2,totfrac3,totfrac4,nchargetot
double precision :: fracijBe,fracijBe2,fracijF,fracijO

double precision,dimension(5,5) :: rdfmin, rdfmin2

double precision, dimension(nmax)   :: x,y,z

character*80 :: filein1
character*2, dimension(5) :: spc

!Reading the input file
open(10,file='speciation.inpt')

read(10,*) nconfigs
read(10,*) nskip
read(10,*) spc(1),nanion1
read(10,*) spc(2),nanion2
read(10,*) spc(3),ncationpol1
read(10,*) spc(4),ncationpol2
read(10,*) spc(5),ncationfree1
read(10,*) boxlength

allocate(spec(0:ncationpol1,0:ncationpol2,0:nanion1,0:nanion2))
allocate(nneighbcat2(ncationpol2))

nion=nanion1+nanion2+ncationpol1+ncationpol2+ncationfree1

rdfmin=0.d0
rdfmin2=0.d0
do i=1,5
   read(10,*)rdfmin(i,1:5)
   do j=1,5
      rdfmin2(i,j)=rdfmin(i,j)**2
   enddo
enddo

read(10,'(a)') filein1

open(11,file=filein1,status='old')
close(10)

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

nchargetot=0
nmol=0
totfrac1=0.0d0
totfrac2=0.0d0
totfrac3=0.0d0
totfrac4=0.0d0

do i=0,ncationpol1
   do j=0,ncationpol2
      do k=0,nanion1
         do l=0,nanion2
            spec(i,j,k,l)=0
         enddo   
      enddo
   enddo
enddo

!loop over nconfigs
do l=1,nskip

!Reading positions.out


do i=1,nion
   read(11,*) x(i),y(i),z(i)
enddo
enddo
do l=1,nconfigs

!Reading positions.out


do i=1,nion
   read(11,*) x(i),y(i),z(i)
enddo




!Initialize variables

do i=1,nmax
   chain(i)=0
   do j=1,nmax
      atoms(i,j)=0
   enddo
enddo

nchain=0
ntri=0

!loop 1 looks for cation1-cation1 linkages

do i=1,ncationpol1

   imark=i+nanion1+nanion2
   nneighb=0
   

   !loop 2 over cation 1
   do j=1,i-1

      jmark=j+nanion1+nanion2

      dx=x(imark)-x(jmark)
      dy=y(imark)-y(jmark)
      dz=z(imark)-z(jmark)

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      !tests  if they are bonded
      !2 tests: d(cation-cation)<rdfmin and presence of 1 or more bridging anion

       if (r2.le.rdfmin2(3,3)) then

          call linkage(imark,jmark,nshar,nanion1,rdfmin2(3,1),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighb=nneighb+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared anions 1 (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
!              if(shar(k).eq.4)write(6,*)imark,jmark
             enddo
          endif

          call linkage2(imark,jmark,nshar,nanion1,nanion2,rdfmin2(3,2),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighb=nneighb+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared anions 2 (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
             enddo
          endif

       endif

    enddo
  
  ! if i has no neighbour then creation of a new chain
    if (nneighb.eq.0) then
       nchain=nchain+1
       chain(imark)=nchain
       atoms(nchain,1)=1
       atoms(nchain,2)=imark
    endif
 enddo

!loop 2 looks for cation1-cation2 linkages

write(6,*)chain(124),'avant boucle 2'
do i=1,ncationpol2

   imark=i+nanion1+nanion2+ncationpol1
   nneighbcat2(i)=0
   

   !loop 2 over cation 1
   do j=1,ncationpol1

      jmark=j+nanion1+nanion2

      dx=x(imark)-x(jmark)
      dy=y(imark)-y(jmark)
      dz=z(imark)-z(jmark)

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      !tests  if they are bonded
      !2 tests: d(cation-cation)<rdfmin and presence of 1 or more bridging anion

       if (r2.le.rdfmin2(4,3)) then

          call linkage3(imark,jmark,nshar,nanion1,rdfmin2(4,1),rdfmin2(3,1),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighbcat2(i)=nneighbcat2(i)+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared anions 1 (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
!              if(shar(k).eq.4)write(6,*)imark,jmark
             enddo
          endif

          call linkage4(imark,jmark,nshar,nanion1,nanion2,rdfmin2(4,2),rdfmin2(3,2),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighbcat2(i)=nneighbcat2(i)+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared anions 2 (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
             enddo
          endif

       endif

    enddo
  
  ! if i has no neighbour then creation of a new chain
    if (nneighbcat2(i).eq.0) then
       nchain=nchain+1
       chain(imark)=nchain
       atoms(nchain,1)=1
       atoms(nchain,2)=imark
    endif
 enddo

write(6,*)chain(124),'avant boucle 3'

!loop 3 looks for cation2-cation2 linkages

do i=1,ncationpol2

   imark=i+nanion1+nanion2+ncationpol1
!  nneighb=0
   

   !loop 2 over cation 2
   do j=1,i-1

      jmark=j+nanion1+nanion2+ncationpol1

      dx=x(imark)-x(jmark)
      dy=y(imark)-y(jmark)
      dz=z(imark)-z(jmark)

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      !tests  if they are bonded
      !2 tests: d(cation-cation)<rdfmin and presence of 1 or more bridging anion

       if (r2.le.rdfmin2(4,4)) then

          call linkage(imark,jmark,nshar,nanion1,rdfmin2(4,1),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighbcat2(i)=nneighbcat2(i)+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared anions 1 (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
!              if(shar(k).eq.4)write(6,*)imark,jmark
             enddo
          endif

          call linkage2(imark,jmark,nshar,nanion1,nanion2,rdfmin2(4,2),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighbcat2(i)=nneighbcat2(i)+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared anions 2 (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
             enddo
          endif

       endif

    enddo
  
  ! if i has no neighbour then creation of a new chain
    if (nneighbcat2(i).eq.0) then
       nchain=nchain+1
       chain(imark)=nchain
       atoms(nchain,1)=1
       atoms(nchain,2)=imark
    endif
 enddo

write(6,*)chain(124),'apres boucle 3'
! Add the unshared anions in the good chains
 do i=1,nanion1

    j=0

    do while (chain(i).eq.0 .and. j.lt.(ncationpol1+ncationpol2))
       
      j=j+1
      jmark=nanion1+nanion2+j
      dx=x(jmark)-x(i) 
      dy=y(jmark)-y(i) 
      dz=z(jmark)-z(i) 

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      if(j.le.ncationpol1)then
         if (r2.le.rdfmin2(3,1)) then
            chain(i)=chain(jmark)
            atoms(chain(i),1)=atoms(chain(i),1)+1
            atoms(chain(i),atoms(chain(i),1)+1)=i
         endif
      else
         if (r2.le.rdfmin2(4,1)) then
            chain(i)=chain(jmark)
            atoms(chain(i),1)=atoms(chain(i),1)+1
            atoms(chain(i),atoms(chain(i),1)+1)=i
         endif
      endif         
    enddo
  ! if F is alone then creates a new chain
     if (chain(i).eq.0) then
        nchain=nchain+1
        chain(i)=nchain
        atoms(nchain,1)=1
        atoms(nchain,2)=i
     endif
 enddo


 ! Add the unshared oxides in the good chains
 do i=nanion1+1,nanion1+nanion2

    j=0

    do while (chain(i).eq.0 .and. j.lt.(ncationpol1+ncationpol2))
       
      j=j+1
      jmark=nanion1+nanion2+j
      dx=x(jmark)-x(i) 
      dy=y(jmark)-y(i) 
      dz=z(jmark)-z(i) 

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      if(j.le.ncationpol1)then
         if (r2.le.rdfmin2(3,2)) then
            chain(i)=chain(jmark)
            atoms(chain(i),1)=atoms(chain(i),1)+1
            atoms(chain(i),atoms(chain(i),1)+1)=i
         endif
      else   
         if (r2.le.rdfmin2(4,2)) then
            chain(i)=chain(jmark)
            atoms(chain(i),1)=atoms(chain(i),1)+1
            atoms(chain(i),atoms(chain(i),1)+1)=i
         endif
      endif   

    enddo
  ! if F is alone then creates a new chain
     if (chain(i).eq.0) then
        nchain=nchain+1
        chain(i)=nchain
        atoms(nchain,1)=1
        atoms(nchain,2)=i
     endif
 enddo






 do i=1,nmax
    nan1(i)=0
    nan2(i)=0
    ncat1(i)=0
    ncat2(i)=0
    ncharge(i)=0
    write(100,*)chain(i)
 enddo

!write(6,*)(atoms(1,i),i=1,atoms(1,1)+1)
! remove triclusters 
 do i=1,nchain
    do j=2,atoms(i,1)
       do k=1,j-1
          if((atoms(i,k+1).ne.0).and.(atoms(i,k+1).eq.atoms(i,j+1)))then
!             write(6,*)atoms(i,k+1),atoms(i,j+1),j,k
!write(6,*)(atoms(1,ll),ll=1,atoms(1,1)+1)
             atoms(i,1)=atoms(i,1)-1
             do ll=k,atoms(i,1)
                atoms(i,ll+1)=atoms(i,ll+2)
             enddo
             atoms(i,atoms(i,1)+2)=0
!             write(6,*)(atoms(1,ll),ll=1,atoms(1,1)+5)
!             atoms(i,atoms(i,1))=0
          endif
       enddo
    enddo
 enddo

 do i=1,nchain
    do j=2,atoms(i,1)
       do k=1,j-1
          if((atoms(i,k+1).ne.0).and.(atoms(i,k+1).eq.atoms(i,j+1)))then
!             write(6,*)atoms(i,k+1),atoms(i,j+1),j,k
!write(6,*)(atoms(1,ll),ll=1,atoms(1,1)+1)
             atoms(i,1)=atoms(i,1)-1
             do ll=k,atoms(i,1)
                atoms(i,ll+1)=atoms(i,ll+2)
             enddo
             atoms(i,atoms(i,1)+2)=0
!             write(6,*)(atoms(1,ll),ll=1,atoms(1,1)+5)
!             atoms(i,atoms(i,1))=0
          endif
       enddo
    enddo
 enddo
!write(6,*)(atoms(1,i),i=1,atoms(1,1)+1)

do i=1,nchain
       do j=1,atoms(i,1)
          if (atoms(i,j+1).le.nanion1) then
             nan1(i)=nan1(i)+1
!            write(6,*)atoms(i,j+1),i
             ncharge(i)=ncharge(i)-2
          elseif((atoms(i,j+1).gt.nanion1).and.(atoms(i,j+1).le.(nanion1+nanion2)))then
             nan2(i)=nan2(i)+1
!             write(26,*)atoms(i,j+1),i,j,atoms(i,1)
             ncharge(i)=ncharge(i)-1
          elseif((atoms(i,j+1).gt.nanion1+nanion2).and.(atoms(i,j+1).le.(nanion1+nanion2+ncationpol1)))then
             ncat1(i)=ncat1(i)+1
             ncharge(i)=ncharge(i)+3
          else    
             ncat2(i)=ncat2(i)+1
             ncharge(i)=ncharge(i)+3
          endif
       enddo
       nchargetot=nchargetot+ncharge(i)
       if ((nan1(i).ne.0).or.(nan2(i).ne.0)) then
          spec(ncat1(i),ncat2(i),nan1(i),nan2(i))=spec(ncat1(i),ncat2(i),nan1(i),nan2(i))+1
          nmol=nmol+1
       endif
 enddo

 do i=1,nchain
    write(101,*)atoms(i,2:atoms(i,1)+1)
 enddo 





 enddo

 open(21,file='speciation-cation1.dat')
 open(22,file='speciation-cation2.dat')
 open(23,file='speciation-anion1.dat')
 open(24,file='speciation-anion2.dat')

 do i=0,ncationpol1
    do j=0,ncationpol2
       do k=0,nanion1
          do l=0,nanion2
             fracijBe=100.0*(float(spec(i,j,k,l)))*float(i)/(float(nconfigs)*float(ncationpol1))
             totfrac1=totfrac1+fracijBe
             if(ncationpol2.ne.0)then
                fracijBe2=100.0*(float(spec(i,j,k,l)))*float(j)/(float(nconfigs)*float(ncationpol2))
                totfrac2=totfrac2+fracijBe2
             endif    
             fracijO=100.0*(float(spec(i,j,k,l)))*float(k)/(float(nconfigs)*float(nanion1))
             totfrac3=totfrac3+fracijO
             if(nanion2.ne.0)then
                fracijF=100.0*(float(spec(i,j,k,l)))*float(l)/(float(nconfigs)*float(nanion2))
                totfrac4=totfrac4+fracijF
             endif   
             if(fracijBe.gt. 0.01) then
               write(21,*) fracijBe,'% de ',spc(3),i,spc(4),j,spc(1),k,spc(2),l
             endif
             if(fracijBe2.gt. 0.01) then
               write(22,*) fracijBe2,'% de ',spc(3),i,spc(4),j,spc(1),k,spc(2),l
             endif
             if(fracijO.gt. 0.01) then
               write(23,*) fracijO,'% de ',spc(3),i,spc(4),j,spc(1),k,spc(2),l
             endif
             if(fracijF.gt. 0.01) then
               write(24,*) fracijF,'% de ',spc(3),i,spc(4),j,spc(1),k,spc(2),l
             endif
          enddo    
       enddo
   enddo
 enddo
 nchargetot=nchargetot/float(nconfigs)
 write(6,*)nchargetot+ncationfree1,totfrac1,totfrac2,totfrac3,totfrac4

 close(21)
 close(22)
 close(23)
 close(24)
 close(11)


 end



subroutine linkage(imark,jmark,nshar,nanion1,rdfcut2,boxlength,x,y,z,shar)

implicit none

integer :: nmax

parameter (nmax=800)

integer :: i,j,k,nshar,nanion1
integer :: imark,jmark
integer, dimension(10) :: shar

double precision, dimension(nmax) :: x,y,z
double precision dx1,dx2,dy1,dy2,dz1,dz2,dr1,dr2,rdfcut2
double precision boxlength,halfbox,halfboxrec

nshar=0

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

do i=1,nanion1

   dx1=x(i)-x(imark)
   dy1=y(i)-y(imark)
   dz1=z(i)-z(imark)

   dx1=dx1-boxlength*int(dx1*halfboxrec)
   dy1=dy1-boxlength*int(dy1*halfboxrec)
   dz1=dz1-boxlength*int(dz1*halfboxrec)

   dx2=x(i)-x(jmark)
   dy2=y(i)-y(jmark)
   dz2=z(i)-z(jmark)

   dx2=dx2-boxlength*int(dx2*halfboxrec)
   dy2=dy2-boxlength*int(dy2*halfboxrec)
   dz2=dz2-boxlength*int(dz2*halfboxrec)

   dr1=dx1**2+dy1**2+dz1**2
   dr2=dx2**2+dy2**2+dz2**2

   if (dr1.le.rdfcut2 .and. dr2.le.rdfcut2) then
      nshar=nshar+1
      shar(nshar)=i
   endif
enddo

return

end

subroutine linkage2(imark,jmark,nshar,nanion1,nanion2,rdfcut2,boxlength,x,y,z,shar)

implicit none

integer :: nmax

parameter (nmax=800)

integer :: i,j,k,nshar,nanion1,nanion2
integer :: imark,jmark
integer, dimension(10) :: shar

double precision, dimension(nmax) :: x,y,z
double precision dx1,dx2,dy1,dy2,dz1,dz2,dr1,dr2,rdfcut2
double precision boxlength,halfbox,halfboxrec

nshar=0

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

do i=nanion1+1,nanion1+nanion2

   dx1=x(i)-x(imark)
   dy1=y(i)-y(imark)
   dz1=z(i)-z(imark)

   dx1=dx1-boxlength*int(dx1*halfboxrec)
   dy1=dy1-boxlength*int(dy1*halfboxrec)
   dz1=dz1-boxlength*int(dz1*halfboxrec)

   dx2=x(i)-x(jmark)
   dy2=y(i)-y(jmark)
   dz2=z(i)-z(jmark)

   dx2=dx2-boxlength*int(dx2*halfboxrec)
   dy2=dy2-boxlength*int(dy2*halfboxrec)
   dz2=dz2-boxlength*int(dz2*halfboxrec)

   dr1=dx1**2+dy1**2+dz1**2
   dr2=dx2**2+dy2**2+dz2**2

   if (dr1.le.rdfcut2 .and. dr2.le.rdfcut2) then
      nshar=nshar+1
      shar(nshar)=i
   endif
enddo

return

end

subroutine linkage3(imark,jmark,nshar,nanion1,rdfcut2i,rdfcut2j,boxlength,x,y,z,shar)

implicit none

integer :: nmax

parameter (nmax=800)

integer :: i,j,k,nshar,nanion1
integer :: imark,jmark
integer, dimension(10) :: shar

double precision, dimension(nmax) :: x,y,z
double precision dx1,dx2,dy1,dy2,dz1,dz2,dr1,dr2,rdfcut2i,rdfcut2j
double precision boxlength,halfbox,halfboxrec

nshar=0

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

do i=1,nanion1

   dx1=x(i)-x(imark)
   dy1=y(i)-y(imark)
   dz1=z(i)-z(imark)

   dx1=dx1-boxlength*int(dx1*halfboxrec)
   dy1=dy1-boxlength*int(dy1*halfboxrec)
   dz1=dz1-boxlength*int(dz1*halfboxrec)

   dx2=x(i)-x(jmark)
   dy2=y(i)-y(jmark)
   dz2=z(i)-z(jmark)

   dx2=dx2-boxlength*int(dx2*halfboxrec)
   dy2=dy2-boxlength*int(dy2*halfboxrec)
   dz2=dz2-boxlength*int(dz2*halfboxrec)

   dr1=dx1**2+dy1**2+dz1**2
   dr2=dx2**2+dy2**2+dz2**2

   if (dr1.le.rdfcut2i .and. dr2.le.rdfcut2j) then
      nshar=nshar+1
      shar(nshar)=i
   endif
enddo

return

end

subroutine linkage4(imark,jmark,nshar,nanion1,nanion2,rdfcut2i,rdfcut2j,boxlength,x,y,z,shar)

implicit none

integer :: nmax

parameter (nmax=800)

integer :: i,j,k,nshar,nanion1,nanion2
integer :: imark,jmark
integer, dimension(10) :: shar

double precision, dimension(nmax) :: x,y,z
double precision dx1,dx2,dy1,dy2,dz1,dz2,dr1,dr2,rdfcut2i,rdfcut2j
double precision boxlength,halfbox,halfboxrec

nshar=0

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

do i=nanion1+1,nanion1+nanion2

   dx1=x(i)-x(imark)
   dy1=y(i)-y(imark)
   dz1=z(i)-z(imark)

   dx1=dx1-boxlength*int(dx1*halfboxrec)
   dy1=dy1-boxlength*int(dy1*halfboxrec)
   dz1=dz1-boxlength*int(dz1*halfboxrec)

   dx2=x(i)-x(jmark)
   dy2=y(i)-y(jmark)
   dz2=z(i)-z(jmark)

   dx2=dx2-boxlength*int(dx2*halfboxrec)
   dy2=dy2-boxlength*int(dy2*halfboxrec)
   dz2=dz2-boxlength*int(dz2*halfboxrec)

   dr1=dx1**2+dy1**2+dz1**2
   dr2=dx2**2+dy2**2+dz2**2

   if (dr1.le.rdfcut2i .and. dr2.le.rdfcut2j) then
      nshar=nshar+1
      shar(nshar)=i
   endif
enddo

return

end
