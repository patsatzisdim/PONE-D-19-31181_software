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
cpd   serum: GLC 5.5, LAC 1  SS values    ! at exersice it changes

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42822d0*V3     
c      y(4)=1.67726d0*V4       
c      y(5)=1.61144d0*V5              
c      y(6)=9.979d-1*V1              
c      y(7)=9.958d-1*V2             
c      y(8)=1.09798d0*V3              
c      y(9)=9.8084d-1*V4              
c      y(10)=8.9903d-1*V5

cpd   serum: GLC 5.5, LAC 8         
c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=3.72262d0*V1              
c      y(7)=1.24981d0*V2             
c      y(8)=1.28519d0*V3              
c      y(9)=1.17320d0*V4              
c      y(10)=1.06285d0*V5
cpd   serum: GLC 5.5, LAC 20                                ! ! ! ! USE THESE ICs FOR EVERY RATE ELEVATED SCENARIO

      y(1)=3.46678d0*V1       
      y(2)=1.96584d0*V2       
      y(3)=1.42821d0*V3     
      y(4)=1.67725d0*V4       
      y(5)=1.61143d0*V5              
      y(6)=6.09233d0*V1              
      y(7)=1.41558d0*V2             
      y(8)=1.41002d0*V3              
      y(9)=1.30144d0*V4              
      y(10)=1.17088d0*V5

cpd   serum: GLC 5.5, LAC 20 + 80% elevated rates R17,R18

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=5.27324d0*V1              
c      y(7)=0.69832d0*V2             
c      y(8)=0.66115d0*V3              
c      y(9)=0.57639d0*V4              
c      y(10)=0.49983d0*V5
c
cpd   serum: GLC 5.5, LAC 20 + 65% elevated rates R17,R18

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=5.35086d0*V1              
c      y(7)=0.76510d0*V2             
c      y(8)=0.73142d0*V3              
c      y(9)=0.64384d0*V4              
c      y(10)=0.56216d0*V5

cpd   serum: GLC 5.5, LAC 10 + 50% elevated rates R17,R18

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=3.79698d0*V1              
c      y(7)=0.77352d0*V2             
c      y(8)=0.77528d0*V3              
c      y(9)=0.68134d0*V4              
c      y(10)=0.60237d0*V5

cpd   serum: GLC 5.5, LAC 20  different paramters (1.5*cr21-cr24) at R^13 

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=6.09261d0*V1              
c      y(7)=1.41583d0*V2             
c      y(8)=1.40993d0*V3              
c      y(9)=1.30148d0*V4              
c      y(10)=1.17091d0*V5

cpd   serum: GLC 5.5, LAC 20  different paramters (1.5*cr21-cr24) at R^13 +80% lactate oxidation rates

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=5.27546d0*V1              
c      y(7)=0.70022d0*V2             
c      y(8)=0.66066d0*V3              
c      y(9)=0.57666d0*V4              
c      y(10)=0.50004d0*V5

cpd   serum: GLC 5.5, LAC 20  different paramters (0.5*cr21-cr24) at R^13 

c      y(1)=3.46678d0*V1       
c      y(2)=1.96584d0*V2       
c      y(3)=1.42821d0*V3     
c      y(4)=1.67725d0*V4       
c      y(5)=1.61143d0*V5              
c      y(6)=6.09165d0*V1              
c      y(7)=1.41498d0*V2             
c      y(8)=1.41023d0*V3              
c      y(9)=1.30136d0*V4              
c      y(10)=1.17082d0*V5

cpd   serum: GLC 5.5, LAC 20  different paramters (1.5*cr5-cr8) at R^3 
c
c      y(1)=3.46744d0*V1       
c      y(2)=1.96684d0*V2       
c      y(3)=1.40688d0*V3     
c      y(4)=1.66781d0*V4       
c      y(5)=1.60219d0*V5              
c      y(6)=6.09156d0*V1              
c      y(7)=1.41490d0*V2             
c      y(8)=1.40926d0*V3              
c      y(9)=1.30076d0*V4              
c      y(10)=1.17030d0*V5

cpd   serum: GLC 5.5, LAC 20  different paramters (0.5*cr5-cr8) at R^3 

c      y(1)=3.46514d0*V1       
c      y(2)=1.96337d0*V2       
c      y(3)=1.48318d0*V3     
c      y(4)=1.70152d0*V4       
c      y(5)=1.63518d0*V5              
c      y(6)=6.09420d0*V1              
c      y(7)=1.41726d0*V2             
c      y(8)=1.41189d0*V3              
c      y(9)=1.30312d0*V4              
c      y(10)=1.17231d0*V5

cpd   serum: GLC 3.33, LAC 10 + 50% elevated rates

c      y(1)=1.94603d0*V1       
c      y(2)=0.86337d0*V2       
c      y(3)=0.41134d0*V3     
c      y(4)=0.59681d0*V4       
c      y(5)=0.55368d0*V5              
c      y(6)=3.75121d0*V1              
c      y(7)=0.72298d0*V2             
c      y(8)=0.71611d0*V3              
c      y(9)=0.63117d0*V4              
c      y(10)=0.55974d0*V5

cpd   For steady-state values set them all zeros
c      do i=1,N
c       y(i)=0.0D0
c      enddo 

      write(*,*) 'Initial Conditions Taken are all zero'
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
  
      return 
      end
cpd------------------------------------------------------------------
