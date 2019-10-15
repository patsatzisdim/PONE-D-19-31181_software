cpd------------------------------------------------------------------
cpd   Initial conditions subroutine
      subroutine InitCond(n,t,y)
      implicit none
      include 'paramet.i'
      integer :: i
      integer, intent(in) :: n
      double precision, intent(out) :: y(n)
      double precision, intent(inout) :: t
      integer :: IOState,iflag


      if (n.lt.10) then
       write(*,*) 'Dimensions of initial conditions changed'
      endif
      iflag=1
      if(iflag.eq.1) then
cpd      
cpd   These are the values of initial conditions for SS
cpd   serum: GLC 5.5, LAC 1  SS values    

      y(1)=3.530786d0*V1       
      y(2)=2.062930d0*V2       
      y(3)=0.9797612d0*V3     
      y(4)=1.757812d0*V4       
      y(5)=1.557273d0*V5              
      y(6)=9.81853d-1*V1              
      y(7)=9.63779d-1*V2             
      y(8)=9.26821d-1*V3              
      y(9)=9.68583d-1*V4              
      y(10)=9.97881d-1*V5

      write(*,*) 'Initial Conditions Taken for steady-state'
      else
cpd   In case of initial conditions from other file
cpd   read from the end of file
      open(3,file='Asol.dat',action='read')
 
      do i=1,100000000
       read(3,*,IOSTAT=IOstate) t,y(1:n)
      if(IOstate.ne.0) exit
      enddo
      close(3)
      write(*,*) 'Initial Conditions Taken from Asol.dat'   
      endif
cpd   For steady-state values set them all zeros
c      do i=1,N
c       y(i)=0.0D0
c      enddo   
      return 
      end
cpd------------------------------------------------------------------
