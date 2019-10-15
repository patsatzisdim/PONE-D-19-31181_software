      subroutine fdjac2(fcn2,n1,n2,x,fvec,fjac2,ldfjac2,epsfcn,wa3)
      integer n1,n2,ldfjac2
C*----*------------------------------------------------------------------*-
c
      double precision epsfcn
      double precision x(n1),fvec(n2),fjac2(ldfjac2,n1),wa3(n2)
ci_start
      external fcn2
ci_stop
c     **********
c
c     subroutine fdjac1
c
c     this subroutine computes a forward-difference approximation
c     to the n by n jacobian matrix associated with a specified
c     problem of n functions in n variables. if the jacobian has
c     a banded form, then function evaluations are saved by only
c     approximating the nonzero terms.
c
c     the subroutine statement is
c
c       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa1)
c                         
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(n,x,fvec,iflag)
c         integer n,iflag
c         double precision x(n),fvec(n)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac1.
c         in this case set iflag to a negative integer.
c
c       n1 is a positive integer input variable set to the number
c         of variables.
c
c       n2 is a positive integer input variable set to the number
c         of functions.
c
c       x is an input array of length n.
c
c       fvec is an input array of length n which must contain the
c         functions evaluated at x.
c
c       fjac is an output n by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than n
c         which specifies the leading dimension of the array fjac.
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac1. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa1 is work array of length n.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar
C*----*------------------------------------------------------------------*-
C                        ... dpmpar replaced by d1mach from blas1 
C                            (Andrew Tron, Princeton University)
C*----*------------------------------------------------------------------*-
c
c       fortran-supplied ... dabs,dmax1,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j
      double precision eps,epsmch,h,temp,zero
C      double precision dpmpar
C*----*------------------------------------------------------------------*-
C     d1mach must be declared as double precision in order to get a
C     valid result  (Andrew Tron, Princeton University)
C*----*------------------------------------------------------------------*-
      double precision d1mach
      data zero /0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = d1mach(4)
c
      eps = dsqrt(dmax1(epsfcn,epsmch))
c
c        computation of dense approximate jacobian.
c
         do 20 j = 1, n1
            temp = x(j)
           h = eps*dmax1(dabs(temp),1.0d0)*dsign(1.0d0,temp)
C
C*----*------------------------------------------------------------------*-
            if (h .eq. zero) h = eps
            x(j) = temp + h
            call fcn2(x,wa3)
C*----*------------------------------------------------------------------*-
C
            x(j) = temp
            do 10 i = 1, n2
               fjac2(i,j) = (wa3(i) - fvec(i))/h
   10          continue
   20       continue
   
c      write(*,*) fjac2(1,1),fjac2(1,2),fjac2(1,3),fjac2(1,4)
c      write(*,*) fjac2(2,1),fjac2(2,2),fjac2(2,3),fjac2(2,4)
c      write(*,*) fjac2(3,1),fjac2(3,2),fjac2(3,3),fjac2(3,4)
c      write(*,*) fjac2(4,1),fjac2(4,2),fjac2(4,3),fjac2(4,4)
c      write(*,*) fjac2(5,1),fjac2(5,2),fjac2(5,3),fjac2(5,4)
c      write(*,*) fjac2(6,1),fjac2(6,2),fjac2(6,3),fjac2(6,4)
 
      return
c
c     last card of subroutine fdjac1.
c
      end
*======================================================================*
*                                              S-STEP comment [15/10/96]
*======================================================================*
*
*  Source:         NETLIB distribution of MINPACK.
*                  http://www.netlib.no/netlib/minpack/fdjac1.f
*                  Modified by Andrew Tron, Princeton University.
*
*======================================================================*

