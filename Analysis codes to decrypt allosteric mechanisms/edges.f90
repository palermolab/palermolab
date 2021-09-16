!==================MAIN PROGRAM===============================!
program edges
implicit none 
REAL,ALLOCATABLE,DIMENSION(:) :: COORD
INTEGER :: natom,nframe,frame,i,j,k,jj,ii,threenatom,INFO,nheavyatom,nres
real :: temp, temp1, temp2, temp3
real,dimension(3) :: temp4,temp5
integer,allocatable,dimension(:) :: INTVEC,VEC,indexheavy,VECSCRATCH
integer,allocatable,dimension(:,:) :: MAT
logical :: ESSENTIAL_DYN,MI 
real :: treshold , rvalue
! MI-damp
logical :: damp
real,allocatable,dimension(:,:) :: distances
! READING VARIABLES
character(2) :: script1
character(1) :: script2
integer, allocatable :: res_index(:),res_index1(:)
character :: algoc1, algoc2
character(5), allocatable ::   resname(:,:)
character(1), allocatable ::   indexname1(:),indexname2(:)
character(5)              ::   algoc3
character(4)              ::   algoc4
character(1)              ::   algo, char5
character(20)             ::   script
integer                   ::   atom_index
real                      ::   treshold_dist
!.............................................................!
!.............................................................!

!*******************! 
!   write(*,*) 'Number of atoms, number of frames '
!   read(*,*) natom,nframe
    write(*,*) 'Reading input'
    open(unit=10,file='edges.inp')
    do i=1,5
       read(10,*) algoc4,rvalue
       if(algoc4.eq.'nres') nres=int(rvalue)
       if(algoc4.eq.'nato') natom=int(rvalue)
       if(algoc4.eq.'nfra') nframe=int(rvalue)
       if(algoc4.eq.'tres') treshold=rvalue
       if(algoc4.eq.'tred') treshold_dist=rvalue
    enddo
    treshold=treshold*nframe
!    nres=1361
!    natom=22436
!    nframe=4224
!    treshold=0.75*nframe
    write(*,*) 'treshold =', int(treshold)
    write(*,*) '# residues =', nres
    write(*,*) '# frames =',nframe
    write(*,*) '# atoms =', natom
    write(*,*) 'treshold distance =', treshold_dist
!   static variables
!*******************!
    allocate(resname(Natom,3))
    allocate(res_index(Natom))
    ALLOCATE(coord(3*natom)) ! The possitions in each frame are stored and in the last possition of coord we store their average

! Read the coordinates in xyz format
    write(*,*) 'Reading coordinates'
    open(1,file='XXX.pdb')

! Primera lectura para obtener nheavyatom    
    ii=1
    jj=0
    read(1,*)
!    read(1,*)
!    read(1,*) 
    do i=1,Natom
        read(1,*) algoc4, atom_index, resname(i,1), resname(i,2), char5, res_index(i), coord(ii), coord(ii+1), coord(ii+2)
!       ATOM      9  HB3 LYS    1      45.275  59.848  51.216  0.00  0.00           H
        ii=ii+3
        write(script1,'(A2)')resname(i,1)
        write(script2,'(A1)')script1
        if(script1.EQ.'1H')res_index(i)=-1
        if(script1.EQ.'2H')res_index(i)=-1
        if(script1.EQ.'3H')res_index(i)=-1
        if(script1.EQ.'4H')res_index(i)=-1
        if(script2.EQ.'H')res_index(i)=-1
        if(res_index(i).ne.-1) jj=jj+1
    enddo
 !   read(1,*)
 !   read(1,*) 
    nheavyatom=jj
    write(*,*) 'number of heavy atoms =', nheavyatom
    ALLOCATE(indexheavy(natom))
    ii=1
    jj=1
    write(*,*) 'calculating the list of heavy atoms'
! now lets find the indexes asociated with the heavy atoms
    do i=1,natom
       if(res_index(i).ne.-1) then
          indexheavy(i)=i
       else
          indexheavy(i)=0
       endif
    enddo
    ii=1
    jj=1
!Now lets start with the real thing
   deallocate(coord)
   deallocate(res_index)
   allocate(res_index(nheavyatom))
   allocate(coord(3*nheavyatom))
   ALLOCATE(INTVEC(nres*(nres+1)/2),VECSCRATCH(nres*(nres+1)/2))
   rewind(1)
   INTVEC=0
   VECSCRATCH=0
!                                                        #  if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!                                                        #  if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
   write(*,*) 'Starting the edges calculation'
   do i=1,nframe
 !   read(1,*)
 !   read(1,*)
    read(1,*)
!    write(*,*) i
        do j=1,natom
           if(indexheavy(j).ne.0) then
             read(1,*) algoc4 , atom_index , resname(j,1), resname(j,2), char5, res_index(jj), coord(ii), coord(ii+1), coord(ii+2)
             ii=ii+3
             jj=jj+1
           else
             read(1,*)
           endif
        enddo
        ii=1
        jj=1
        do j=1,nheavyatom
           do k=j,nheavyatom
              temp=dist(coord,j,k)
              if((temp.le.treshold_dist).and.(res_index(k).ne.(res_index(j)+1)).and.(res_index(k).ne.(res_index(j)-1)).and.(res_index(k).ne.res_index(j)) ) then
                 VECSCRATCH(res_index(j) + (res_index(k)-1)*res_index(k)/2)=1
              endif
           enddo
        enddo
        INTVEC=INTVEC+VECSCRATCH
        VECSCRATCH=0
 !       read(1,*)
 !       read(1,*)
   enddo
   DEALLOCATE(VECSCRATCH)

   write(*,*) 'evaluation of exclusion list'
!! now lets evaluate the exclusion list for each node
  open(unit=22,file='notfour.txt')
  temp=0
  jj=1
  ii=0
!  do i=1,nres
!     do j=i+1,nres
!        if(INTVEC(i + (j-1)*j/2).lt.int(treshold)) ii=ii+2
!     enddo
!     ii=ii+1  ! for the diagonal
!     ALLOCATE(VEC(ii))
!     do j=1,nres
!        if((j.lt.i).and.(INTVEC(i + (j-1)*(2*nres-j)/2).lt.int(treshold))) then
!           VEC(jj)=j
!           jj=jj+1
!        endif
!        if ((j.ge.i).and.(INTVEC(i + (j-1)*j/2).lt.int(treshold))) then
!          VEC(jj)=j
!          jj=jj+1
!        endif
!     enddo
!     jj=1
!     write(22, '(1000I4)') (vec(j) ,j=1,ii)
!     ii=0
!     DEALLOCATE(VEC)
!  enddo
!
  ALLOCATE(MAT(nres,nres))
  call spunpack(nres,INTVEC,MAT)
  open(unit=20,file='distance-mat.txt')
  do i=1,nres
     write(20,*) ''
     do j=1,nres
        if((j.eq.i).or.(j.eq.i-1).or.(j.eq.i+1)) then
           write(20,*) real(i),real(j), real(nframe)
        else
           write(20,*) real(i),real(j),real(MAT(i,j))
        endif
     enddo
  enddo
  open(unit=22,file='notfour.txt')
  temp=0
  jj=1
  ii=0
  do i=1,nres
     do j=1,nres
        if(MAT(i,j).lt.int(treshold)) ii=ii+1
     enddo
     ALLOCATE(VEC(ii))
     do j=1,nres
        if(MAT(i,j).lt. int(treshold)) then
           VEC(jj)=j
           jj=jj+1
        endif
     enddo
     jj=1
     write(22, '(2000I5)') (vec(j) ,j=1,ii)
!     write(22,*) (vec(j) ,j=1,ii)
     ii=0
     DEALLOCATE(VEC)
  enddo
! 
  contains
!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
       REAL function  DIST(coord,atom1,atom2) 
       implicit none
       integer,intent(in)    :: atom1,atom2
       real, intent(in)      :: coord(:)
       integer               :: ii,jj
       real                  :: temp(3)
       
       ii=(atom1-1)*3+1
       jj=(atom2-1)*3+1
       temp(1) = coord(ii)-coord(jj)
       temp(2) = coord(ii+1)-coord(jj+1)
       temp(3) = coord(ii+2)-coord(jj+2)
       DIST=NORM2(temp)

       end function DIST
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
       SUBROUTINE spunpack(NM,Vector,Matrix)
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: NM
       INTEGER,INTENT(IN)  :: Vector(NM*(NM+1)/2)
       INTEGER,INTENT(OUT) :: Matrix(NM,NM)
       INTEGER            :: ii,jj,idx
 
         DO jj=1,NM
         DO ii=1,jj
           idx=ii+(jj*(jj-1)/2)
           Matrix(ii,jj)=Vector(idx)
           Matrix(jj,ii)=Vector(idx)
         ENDDO
         ENDDO


       RETURN;END SUBROUTINE
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  end program
