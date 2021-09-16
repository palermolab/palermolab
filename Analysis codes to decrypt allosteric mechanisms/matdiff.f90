!==================MAIN PROGRAM===============================!
program matdiff
implicit none 
real :: algo1,algo2,scratch1,scratch2,scratch3,scratch4,sd
integer :: i,j,k,nres
real, allocatable :: vector(:), vector1(:), vector2(:), vector3(:)
  nres=121
  allocate(vector(nres),vector1(nres),vector2(nres),vector3(3*nres))
  vector=0
  vector1=0
  vector2=0
  open(unit=18,file='wt-dist-mat.txt')
  open(unit=19,file='855-dist-mat.txt')
  open(unit=20,file='848-dist-mat.txt')
  open(unit=21,file='810-dist-mat.txt')
  open(unit=22,file='855-wt-diff-mat.txt')
  open(unit=23,file='848-wt-diff-mat.txt')
  open(unit=24,file='810-wt-diff-mat.txt')
  do i=1,nres
     read(18,*)
     read(19,*)
     read(20,*) 
     read(21,*)
     write(22,*) ''
     write(23,*) ''
     write(24,*) ''
     do j=1,nres
        if((j.eq.i).or.(j.eq.i-1).or.(j.eq.i+1)) then
           read(18,*) algo1, algo2, scratch1
           read(19,*) algo1, algo2, scratch2
           read(20,*) algo1, algo2, scratch3
           read(21,*) algo1, algo2, scratch4
           write(22,*) real(i),real(j), scratch2-scratch1
           write(23,*) real(i),real(j), scratch3-scratch1
           write(24,*) real(i),real(j), scratch4-scratch1
           vector(i)=vector(i)+abs(scratch2-scratch1)
           vector1(i)=vector1(i)+abs(scratch3-scratch1)
           vector2(i)=vector2(i)+abs(scratch4-scratch1)
        else
           read(18,*) algo1, algo2, scratch1
           read(19,*) algo1, algo2, scratch2
           read(20,*) algo1, algo2, scratch3
           read(21,*) algo1, algo2, scratch4
           write(22,*) real(i),real(j), scratch2-scratch1
           write(23,*) real(i),real(j), scratch3-scratch1
           write(24,*) real(i),real(j), scratch4-scratch1
           vector(i)=vector(i)+abs(scratch2-scratch1)
           vector1(i)=vector1(i)+abs(scratch3-scratch1)
           vector2(i)=vector2(i)+abs(scratch4-scratch1)
        endif
     enddo
  enddo
  open(unit=25,file='855-WT-contact-distorsion.txt')
  open(unit=26,file='848-WT-contact-distorsion.txt')
  open(unit=27,file='810-WT-contact-distorsion.txt')
  scratch1=0
  scratch2=0
  scratch3=0
  DO i=1,nres
     vector3(i)=vector(i) 
     vector3(nres+i)=vector1(i)
     vector3(2*nres+i)=vector2(i)
  ENDDO
  call avg(vector3,3*nres,scratch4,sd)
  scratch4=scratch4+sd
  do i=1,nres
     write(25,*) vector(i)
     write(26,*) vector1(i)
     write(27,*) vector2(i)
     if(vector(i).ge.scratch4) scratch1=scratch1+vector(i)
     if(vector1(i).ge.scratch4) scratch2=scratch2+vector1(i)
     if(vector2(i).ge.scratch4) scratch3=scratch3+vector2(i)
  enddo

  open(unit=30,file='accumulated-distorsion.txt')
  write(30,*) '855  848  810'
  write(30,*) scratch1/nres, scratch2/nres,scratch3/nres
  contains
  subroutine avg(vec,N,average,sd)
  implicit none
  integer :: N
  real :: vec(N), average,sd
  
  sd=0.00000
  average=sum(vec)
  average=average/N
  do i=1,N
     sd=sd+(vec(i)-average)**2
  ENDDO
  sd=sd/N
  sd=sqrt(sd)
  write(*,*) 'average =', average
  write(*,*) 'standard deviation =', sd
  end subroutine
  end program
