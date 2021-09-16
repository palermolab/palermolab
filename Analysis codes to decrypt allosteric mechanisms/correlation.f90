!==================MAIN PROGRAM===============================!
program linear_correlation
implicit none 
REAL,ALLOCATABLE,DIMENSION(:,:) :: COORD
REAL,ALLOCATABLE,DIMENSION(:,:) :: COVARMAT,TWOBODYCORR,EIGENVEC
REAL,ALLOCATABLE,DIMENSION(:)   :: COVARVEC,EIGENVAL,WORK
INTEGER :: natom,nframe,frame,i,j,k,jj,ii,threenatom,INFO
real :: temp, temp1, temp2, temp3
logical :: ESSENTIAL_DYN,MI  
! Three body variables
logical :: threebody
REAL,ALLOCATABLE,DIMENSION(:) :: MARGINAL
REAL,ALLOCATABLE,DIMENSION(:,:) :: AVERAGEPROD 
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: COVARMAT3B,THREEBODYCORR
!Centrality
logical :: centrality
REAL,ALLOCATABLE,DIMENSION(:) ::  TWOBODYVEC
! RMSF
logical :: dormsf
real    :: rmsf 
real    :: rmtemp,rmtemp0
! MI-damp
logical :: damp
real,allocatable,dimension(:,:) :: distances
! READING VARIABLES
character (4)     :: char1
integer           :: int2
character (2)     :: char3
character (3)     :: char4
integer           :: int5
!.............................................................!
! linear correlation program calculates:
! (1) The essential dynamics of a trajectory, i.e. the covariance matrix and its eigenvalues and eigenvectors (covar_matrix.txt, covar_eigenval.txt, covar_eigenvec.txt)
! (2) The linear generalized correlation coefficient matrix (mutual information matrix, see [Lange, O. F. and Grubmüller, H. (2006). Proteins, 62: 1053–1061]) and its eigenvalues and eigenvectors.
! (3) Eigenvector centrality coefficients can be calculated for both the essential dynamics and the mutual information matrix.
! (4) Three body correlation matrix.
! Important notes: 
! . Needs to be linked to LAPACK.
!   To do this, you may need to set:
!   LDFLAGS:  -L/usr/local/opt/lapack/lib
!   CPPFLAGS: -I/usr/local/opt/lapack/include
! . It is important to center the trajectory before running this code. That can be done with VMD.
! . This code is not optimized at all!!!! This is just a toy code that will be parallelized and optimized in the future.   
!.............................................................!

! Harcoded logicals !
   MI=.true.           ! calculates mutual information (2 body correlation) matrix
   ESSENTIAL_DYN=.true. ! calculates the eigenvalues/eigenvectors of the covariance matrix
   centrality=.true.
   threebody=.false.
   dormsf=.true.
   damp=.false.   ! long-short range filter 
!*******************! 
!   write(*,*) 'Number of atoms, number of frames '
!   read(*,*) natom,nframe
    natom=1   ! number of atoms in the system
    nframe=1  ! number of frames in the trajectory
! static variables
   threenatom=3*natom
!*******************!
   ALLOCATE(coord(3*natom,nframe+1),COVARMAT(3*natom,3*natom),TWOBODYCORR(natom,natom)) ! The possitions in each frame are stored and in the last possition of coord we store their average
                                                                                        ! COVARMAT and TWOBODYCORR are the covariance matrix and the mutual information matrix respectively

! Read the coordinates in xyz format
   write(*,*) 'Reading coordinates'
   open(1,file='NAME.pdb') ! NAME OF TRAJ FILE
   do frame=1,nframe
      ii=1
      do i=1,natom
! read(1,*) char1,int2,char3 , char4, int5, coord(ii,frame), coord(ii+1,frame), coord(ii+2,frame)
         read(1,*) char1,int2,char3 , char4, int5, coord(ii,frame), coord(ii+1,frame), coord(ii+2,frame)
         ii=ii+3
!!!!!!!!!!         ATOM  304   CA      ILE    304       4.243          23.573             -5.248               0.00           0.00           C
      enddo
   enddo

! Calculate the average coordinate 
  write(*,*) 'Calculating averages'
  call average(coord,natom,nframe) ! the average matrix is stored in the coord(:,nframe+1)

  if(damp) then
    allocate(distances(natom,natom))
    write(*,*) 'Calculating average distances'
    jj=1
    ii=1
    do i=1,natom
       do j=1,natom
          temp=(coord(ii,nframe+1)-coord(jj,nframe+1))**2 + (coord(ii+1,nframe+1)-coord(jj+1,nframe+1))**2 + (coord(ii+2,nframe+1)-coord(jj+2,nframe+1))**2
          temp=sqrt(temp)
          jj=jj+3
          distances(i,j)=temp
          write(50,*) temp
       enddo
       ii=ii+3
       jj=1
    enddo
    jj=1
    ii=1
  endif
  if(DORMSF) then
     rmsf=0.0000000
     write(*,*) 'Calculating RMSF'
     open(33,file='./RMSF.dat')
     rmsf=0
     ii=1
     do i=1,natom
        temp= coord(ii,nframe+1)
        temp1= coord(ii+1,nframe+1)
        temp2= coord(ii+2,nframe+1)
        rmsf=0.000000
        do j=1,nframe
           rmtemp = (coord(ii,j)-temp)**2 + (coord(ii+1,j)-temp1)**2 + (coord(ii+2,j)-temp2)**2
           rmsf=rmsf+rmtemp/(nframe)
        enddo
        ii=ii+3
        write(33,*) sqrt(rmsf)
        rmsf=0.00000000
     enddo
     temp=0.00000
     temp1=0.0000
     temp2=0.0000
  endif

!! The coordenate vector is transformed into a displacement vector substracting the average from the coordinates i.e.  X = X - <X>
   write(*,*) 'Evaluating displacement vector (X - <X>)'
   do i=1,3*natom
      do j=1,nframe
         coord(i,j)=coord(i,j)-coord(i,nframe+1)
      enddo
   enddo 

!! Now lets calculate the covariance matrix

   write(*,*) 'Evaluating covariance matrix'   
   COVARMAT=COVAR(COORD,natom,nframe)         ! COVAR function is the primary target to be optimized!!!!!!!

!! Now lets calculate the essential modes
   
   if(ESSENTIAL_DYN) then
      allocate(COVARVEC(threenatom*(threenatom+1)/2),EIGENVAL(threenatom),EIGENVEC(threenatom,threenatom),WORK(3*threenatom))
      do i=1,3*natom
           do j=i,3*natom
              COVARVEC(i + (j-1)*j/2) = COVARMAT(i,j)     ! This is probably not needed...
           enddo
      enddo
      call SSPEV('V','U',threenatom,COVARVEC,EIGENVAL,EIGENVEC,threenatom,WORK,INFO)
!........................................!
!SSPEV computes all the eigenvalues and, optionally, eigenvectors of a real symmetric matrix A in packed storage.
!subroutine dspev	(	character 	JOBZ,    # JOBZ is CHARACTER*1 = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character 	UPLO,                                    # UPLO is CHARACTER*1 = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
!integer 	N,                                       # N is INTEGER. The order of the matrix A.  N >= 0.
!single  precision, dimension( * ) 	AP,              # AP is SINGLE PRECISION array, dimension (N*(N+1)/2)
!                                                        # On entry, the upper or lower triangle of the symmetric matrix
!                                                        # A, packed columnwise in a linear array.  The j-th column of A
!                                                        # is stored in the array AP as follows:
!                                                        # if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!                                                        # if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!                                                        # On exit, AP is overwritten by values generated during the
!                                                        # reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!                                                        # and first superdiagonal of the tridiagonal matrix T overwrite
!                                                        # the corresponding elements of A, and if UPLO = 'L', the
!                                                        # diagonal and first subdiagonal of T overwrite the
!                                                        # corresponding elements of A.
!single precision, dimension( * ) 	W,               # W is REAL array, dimension (N) If INFO = 0, the eigenvalues in ascending order.
!single precision, dimension( ldz, * ) 	Z,               # Z is REAL array, dimension (LDZ, N)
!                                                        # If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!                                                        # eigenvectors of the matrix A, with the i-th column of  holding the eigenvector associated with W(i).
!                                                        # If JOBZ = 'N', then Z is not referenced.
!integer 	LDZ,                                     # LDZ (INTEGER) is the leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
!single precision, dimension( * ) 	WORK,            # WORK is REAL array, dimension (3*N)
!integer 	INFO                                     # INFO is INTEGER
!                                                        # = 0:  successful exit.
!                                                        # < 0:  if INFO = -i, the i-th argument had an illegal value.
!                                                        # > 0:  if INFO = i, the algorithm failed to converge; i
!                                                        # off-diagonal elements of an intermediate tridiagonal
!                                                        # form did not converge to zero.
!)	
!........................................!
! check if the diagonalization went ok
      if(INFO.lt.0) then
        write(*,*) 'the', INFO ,'-th argument had an illegal value.'
        stop
      endif
      if(INFO.gt.0) stop 'The covariance matrix could not be diagonalized'

! Write the essential modes 
      open(2,file='covar_eigenvec.txt')
      open(3,file='covar_eigenval.txt')
      open(4,file='covar_norm_centrality.txt')
      do i=1,threenatom
            write(2,*) '<==============>'
            write(3,*) i,EIGENVAL(i)
         do j=1,threenatom
            write(2,*) i,EIGENVEC(j,i)
         enddo
      enddo
      ii=1
         write(4,*) EIGENVAL(3*natom)
      do i=1,natom
         write(4,*) sqrt(EIGENVEC(ii,3*natom)**2+EIGENVEC(ii+1,3*natom)**2+EIGENVEC(ii+2,3*natom)**2)
         ii=ii+3
      enddo
      DEALLOCATE(COVARVEC,EIGENVEC,EIGENVAL,WORK) ! Free some memory space
   endif
   

   if(MI) then   
!! Now lets calculate the Mutual information (2 body correlation) matrix
      write(*,*) 'Evaluating mutual information matrix'
      if(threebody) ALLOCATE(MARGINAL(natom))
      TWOBODYCORR=0.000000000
      do i=1,natom
         temp=DETv1(COVARMAT,i,i,3)
         if(threebody) MARGINAL(i)=temp
         do j=1,i-1
            temp1=DETv1(COVARMAT,j,j,3)
            temp2=DETv1(COVARMAT, i,j,6) 
            temp3=temp*temp1/temp2
            temp3=LOG(temp3)
            TWOBODYCORR(i,j)=temp3
            TWOBODYCORR(j,i)=TWOBODYCORR(i,j)
         enddo
      enddo 
      TWOBODYCORR=0.500000*TWOBODYCORR
   endif

   if(MI) then
!! Now lets transform the 2 body correlation into the generalized correlation
!coefficients defined as:
!         c(ij)=[1-exp(-2I{i,j}/3)],
! where I{i,j} is the mutual information or two body correlation we have just
! determined.
      do j=1,natom
         do i=1,j-1
            TWOBODYCORR(i,j)=exp((-0.666666666666667)*TWOBODYCORR(i,j))
            TWOBODYCORR(i,j)=1.000000-TWOBODYCORR(i,j)
            TWOBODYCORR(i,j)=sqrt(TWOBODYCORR(i,j))
            TWOBODYCORR(j,i)=TWOBODYCORR(i,j)
         enddo
         TWOBODYCORR(j,j)=1.00000000000
      enddo
      open(5,file='generalized_correlation_coeff.txt')
      do i=1,natom
         write(5,*)''
         do j=1,natom
            write(5,*) real(j,8),real(i,8),TWOBODYCORR(j,i)
         enddo
      enddo
      open(unit=787,file='dist.txt')   ! WRITING CORRELATION MATRIX (SEE PNAS RIVALTA ET AL 2012)
      write(787,*) natom
      do i=1,natom
        write(787,'(2000F10.6)') (TWOBODYCORR(i,j),j=1,natom)
      enddo
      if(centrality) then
         ALLOCATE(TWOBODYVEC(natom*(natom+1)/2),EIGENVEC(natom,natom),EIGENVAL(natom),WORK(3*natom))
         do i=1,natom
             TWOBODYVEC(i + (i-1)*i/2)=0.0000000    ! Metodo de Christian  
             do j=i+1,natom
                TWOBODYVEC(i + (j-1)*j/2) = TWOBODYCORR(i,j)     ! This is probably not needed...
             enddo
         enddo
         call SSPEV('V','U',natom,TWOBODYVEC,EIGENVAL,EIGENVEC,natom,WORK,INFO)
! check if the diagonalization went ok
      if(INFO.lt.0) then
        write(*,*) 'the', INFO ,'-th argument had an illegal value.'
        stop
      endif

      if(INFO.gt.0) stop 'The covariance matrix could not be diagonalized'
         open(12,file='centrality_eigenvec.txt')    ! THIS IS THE CENTRALITY
         open(13,file='centrality_eigenval.txt')    ! 
!         write(12,*) EIGENVAL(natom)
         do i=1,natom
             write(13,*) i,EIGENVAL(i)
             write(12,*) EIGENVEC(i,natom)
         enddo
         DEALLOCATE(TWOBODYVEC,EIGENVEC,EIGENVAL,WORK) ! Free some memory space
      endif
   endif
contains
!====================AUXILIAR FUNCTIONS======================!
      function COVAR(COORD,natom,nframe) result(COVARMAT)
      implicit none
!....................................................................................................!
!     Covariance matrix calculation:
!     This is maybe the simplest and most intuitive way to calculate the covariance matrix.   
!     This is currently the most demanding computation in the code and can be improved a lot.
!....................................................................................................!
      real, intent(in)    ::   COORD(:,:)
      integer, intent(in) ::   natom
      integer, intent(in) ::   nframe
      real,allocatable    ::   COVARMAT(:,:)
      real                ::   scratch,oneovernframe
      integer             ::   N,k,i,j 
     
      allocate(COVARMAT(3*natom,3*natom))
      oneovernframe=1.00000/nframe 
      scratch=0.00000
      do i=1,3*natom
         do j=1,i
           do k=1,nframe
               scratch=scratch+coord(i,k)*coord(j,k)
           enddo
           scratch=scratch*oneovernframe
           COVARMAT(i,j)=scratch
           COVARMAT(j,i)=scratch
           scratch=0.00000000
         enddo
       enddo
       return;end function 
!!------------------------------------------------------------!
       subroutine  AVERAGE(coord,natom,nframe) 
       implicit none
!....................................................................................................!
!      This subroutine calculates the average value of the coordinates and stores those values in coord(:,nframe+1)
!....................................................................................................!
       integer,intent(in)    :: natom,nframe
       real, intent(inout)      :: coord(3*natom,nframe+1)
       real                  :: avg,oneovernframe
       integer               :: i,j
       
       oneovernframe=1.0000/nframe
       do i=1,3*natom
          coord(i,nframe+1)=0.00000 
          avg=sum(coord(i,:))
          avg=avg/nframe
          coord(i,nframe+1)=avg
       enddo
     
       return;end subroutine
!!------------------------------------------------------------!
      function DET(COVARMAT,i,j) result(determinant)
       implicit none
!....................................................................................................!
! 3x3 determinant:
!      |a  b  c| 
!  det |d  e  f| = aei + bfg + cdh - ceg - bdi - afh 
!      |g  h  i|
! Carefull must be taken, i and j are indexes that label each atom while the indexes in COVARMAT span over al the coordenates of each atom.
!....................................................................................................!
       integer, intent(in)   ::  i,j
       real, intent(in)      ::  COVARMAT(:,:)
       real                  ::  determinant
       real                  ::  temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+3)
       determinant=temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+1)
       determinant=determinant+temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+2)
       determinant=determinant+temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+1)
       determinant=determinant-temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+3)
       determinant=determinant-temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+2)
       determinant=determinant-temp
      
      end function DET
!------------------------------------------------------------!
    REAL function DETv1(COVARMAT, ii,jj,n)
    IMPLICIT NONE
!.....................................................................................................!
! ONLY WORKS WHEN n=6 or n=3 (i.e. 6x6 or 3x3 matrices) !!!!!!!!!!!!!
! 3x3 or 6x6 determinant
! The subroutine is based on two key points:
! [1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!     row operations (column operations would work as well) are used
!     to convert the matrix into upper traingular form
! [2] The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
! Note: In order to make this function more efficient, it would be better to get rid of the variable matrix or vectorize the loops.
!.....................................................................................................!
    INTEGER, INTENT(IN)  :: n,ii,jj
    REAL, INTENT(IN)     :: COVARMAT(:,:)
    REAL, DIMENSION(n,n) :: matrix
    REAL                 :: m, temp
    INTEGER              :: i, j, k, l
    LOGICAL              :: DetExists = .TRUE.
    l = 1
    if(n.eq.6) then
       do i=1,3
          do j=1,3 
             matrix(i,j)=COVARMAT((ii-1)*3+i,(ii-1)*3+j)
             matrix(i,j+3)=COVARMAT((ii-1)*3+i,(jj-1)*3+j)
             matrix(i+3,j)=COVARMAT((jj-1)*3+i,(ii-1)*3+j)
             matrix(i+3,j+3)=COVARMAT((jj-1)*3+i,(jj-1)*3+j)
          enddo
       enddo
    else if(n.eq.3) then
       do i=1,3
          do j=1,3
             matrix(i,j)=COVARMAT((ii-1)*3+i,(jj-1)*3+j)
          enddo
       enddo
    else 
       stop 'wrong determinant size, the current version of DETv1 only allows determinant of 6x6 or 3x3 matrices'
    endif
    
    !Convert to upper triangular form
    do k = 1, n-1
        if (matrix(k,k) == 0) then
            DetExists = .FALSE.
            do i = k+1, n
                if (matrix(i,k) /= 0) then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    enddo
                    DetExists = .TRUE.
                    l=-l
                    exit
                endif
            end do
            if (DetExists .EQV. .FALSE.) then
               DETv1=0.0000
               write(15,*) n,ii,jj
               return
            endif
        endif
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            enddo
        enddo
    enddo
    
    !Calculate determinant by finding product of diagonal elements
    DETv1 = l
    do i = 1, n
        DETv1 = DETv1 * matrix(i,i)
    enddo
    
end function DETv1

    REAL function DETv2(COVAR,n)
    IMPLICIT NONE
!.....................................................................................................!
! ONLY WORKS WHEN n=6 or n=3 (i.e. 6x6 or 3x3 matrices) !!!!!!!!!!!!!
! 3x3 or 6x6 determinant
! The subroutine is based on two key points:
! [1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!     row operations (column operations would work as well) are used
!     to convert the matrix into upper traingular form
! [2] The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
! Note: In order to make this function more efficient, it would be better to get rid of the variable matrix or vectorize the loops.
!.....................................................................................................!
    INTEGER, INTENT(IN)  :: n
    REAL,intent(in)      :: COVAR(:,:)
    REAL,allocatable     :: matrix(:,:)
    REAL                 :: m, temp
    INTEGER              :: i, j, k, l
    LOGICAL              :: DetExists = .TRUE.
    l = 1
    if(n.eq.6) then
       allocate(matrix(n,n))
       matrix=covar
    elseif(n.eq.3) then
       allocate(matrix(n,n))
       do i=1,3
          do j=1,3
             matrix(i,j)=COVAR(i+3,j+3)
          enddo
       enddo
    else 
       stop 'wrong determnant dimension'
    endif
    !Convert to upper triangular form
    do k = 1, n-1
        if (matrix(k,k) == 0) then
            DetExists = .FALSE.
            do i = k+1, n
                if (matrix(i,k) /= 0) then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    enddo
                    DetExists = .TRUE.
                    l=-l
                    exit
                endif
            end do
            if (DetExists .EQV. .FALSE.) then
               DETv2=0.0000
               stop 'Determinants could not be calculated'
            endif
        endif
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            enddo
        enddo
    enddo

    !Calculate determinant by finding product of diagonal elements
    DETv2 = l
    do i = 1, n
        DETv2 = DETv2 * matrix(i,i)
    enddo

    end function DETv2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    end program
