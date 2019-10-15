      program Sol_comparison
      implicit none

      integer, parameter :: N=10
      double precision :: T, asolm1(N), asolm2(N), dum1, dum2, rel(N)
      integer :: i


      open(11,file='Asolm.dat',form='formatted',action='read')
      open(12,file='Asolm2.dat',form='formatted',action='read')
      open(13,file='RelAsol.dat',form='formatted',status='unknown')

      T=0.0d0
      do while(T.lt.3.50d+3) 
       read(11,30) T, (asolm1(i),i=1,N), dum1, dum2
       read(12,30) T, (asolm2(i),i=1,N), dum1, dum2
       do i=1,N
       	rel(i)=dabs(asolm1(i)-asolm2(i))/asolm1(i)
       enddo
       write(13,31) T, (rel(i),i=1,N)
      enddo

   30 format(E18.11,12(E25.16))
   31 format(E18.11,10(E25.16))
      end