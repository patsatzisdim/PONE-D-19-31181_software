c   
c     This sub is to eventually change the paramet.i file
c     by changing the parameters 
c
c     Inputs:
c     iln:     the line of the paramet to be changed in paramet.i file
c     varName:     the name of variable at this line (it is not nessecery, but i will fix it later on)   
      subroutine ChangeParamet(iln,varName,newVal)
      implicit none
      integer, intent(in) ::iln
      double precision, intent(in) :: newVal
      character (len=*) :: varName
      double precision :: varVal
      character :: ch1*17,ch2*1,ch3*1
      integer :: i

      ch1='      parameter ('
      ch2='='
      ch3=')'
      open(1,file='paramet.i',access='sequential',
     - action='readwrite',form='formatted')
      do i=1,100 
       if (i.eq.iln) then 
        write(1,10) ch1,varName,ch2,newVal,ch3
       else
        read(1,*)
       endif
      enddo
c      write(*,10) ch1,varName,ch2,newVal,ch3
   10 format(a17,a3,a1,D10.4,a1)


      return
      end