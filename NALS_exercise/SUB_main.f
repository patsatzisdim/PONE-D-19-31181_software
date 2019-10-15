      program main
      implicit none

      EXTERNAL FEX, JEX 
      include 'paramet.i'

      character sp      
      integer,parameter :: n=10 
      integer,parameter :: m=30     
      integer :: i,k,j
c      INTEGER :: I,J,INFO,ITEMP      

cpd---------------------------------------LSODE work parameters 
      integer :: MF,ITOL,ISTATE,ITASK,IOPT,LRW,LIW
      double precision :: RTOL,ATOL
      double precision, allocatable, dimension(:):: RWORK
      integer, allocatable, dimension(:):: IWORK
cpd---------------------------------------
      double precision, dimension(n) :: yy,ydot,rlmod
      double precision, dimension(n+2) :: yym ! molar concentration + total concentrations
      double precision :: t1,t,dtStart,dt,tout,tend,tend2      
      double precision, dimension(n,n) :: pd,rjac,alpha,betta,test,po,
     - alphaCSP,bettaCSP,SlowPO
      double precision, dimension(:,:), allocatable :: Ar1,Br1,As1,Bs1,
     - Ar2,Br2,As2,Bs2
      integer :: ML, MU, NRPD
      double precision, dimension(n,2) :: reim,tau
      double precision :: tpi(n,m),api(n,m),fip(n),RR(m),impi(n,m)
      character :: varName*3
      double precision :: CSPrtol,CSPatol,temp,sum1
      integer :: noEM,index
      double precision :: PoS(n,n), tpiS(n,m), apiS(n,m), impiS(n,m)
      integer :: PoI(n,n), tpiI(n,m), apiI(n,m), impiI(n,m)
      double precision :: Ras,Rne
      common k
      k=m ! the number of Reactions Rates
            
c      DOUBLE PRECISION :: T,DT,DTstart=1.0d-5  TOUT,
	  ! VECTOR OF VARIABLES : Y(N) 
      
c      DOUBLE PRECISION, DIMENSION(N) :: F1
          ! SOLVING PARAMETERS     
c      DOUBLE PRECISION, DIMENSION(N,N) :: PD,DP,RJAC,ALPHA,BETTA,REIM
c      DOUBLE PRECISION, DIMENSION(N) :: PN,RLMOD
c      DOUBLE PRECISION :: qu1,qu2,EPS
c      DOUBLE PRECISION, DIMENSION(N,M) :: AR,AR1,AS,AS1,AS2	
c      DOUBLE PRECISION, DIMENSION(M,N) :: BR,BR1,BR2,BS,BS1
c      INTEGER, DIMENSION(N) :: IPVT
          ! APROXIMATION VARIABLES
c      DOUBLE PRECISION, DIMENSION(N) :: SMX,SMY,SMPEA
c      DOUBLE PRECISION :: MX,MY,MPEA,Del1,Del2

cpd--------------------------------------------------
      sp=char(9) ! necessary gap used for write	  		  
cpd--------------------------------------------------  
      t=0.0d0 
cpd--------------------------------------------------
cpd   connect subroutine to get the initial conditions
      call InitCond(n,t,yy)    
cpd--------------------------------------------------
cpd   time and more
      dtStart=1.0d-5   
      dt=dtStart
      tout=t+dt
      tend=t+1.02d+3
      write(*,*) "Running till time t=",tend
cpd--------------------------------------------------      
cpd   Lsode parameters, check lsode documentation. 
      MF=21	  
      LRW=22+9*n+n**2
      LIW=20+n
      ALLOCATE(RWORK(LRW))
      ALLOCATE(IWORK(LIW))
      ITOL=1
      RTOL=5.0D-16
      ATOL=1.0D-18
      ITASK=1
      ISTATE=1
      IOPT=1                      ! optional inputs (0.0d0 or 0 means the default)
      RWORK(5)=0.0D0              ! initial stepsize H0
      RWORK(6)=0.0D0              ! maximum stepsize allowed
      RWORK(7)=0.0D0              ! minimum stepsize allowed
      RWORK(8)=0.0D0              !
      RWORK(9)=0.0D0              !
      RWORK(10)=0.0D0             ! 
      IWORK(5)=0                  ! maximum order (integration method)
      IWORK(6)=50000000           ! maximum number of internal steps 
      IWORK(7)=100                ! maximum number of problems printed
      IWORK(8)=0                  !
      IWORK(9)=0                  !
      IWORK(10)=0                 !   
cpd--------------------------------------------------
cpd   open .dat files to write results
      open(101,file='Asol.dat',form='formatted',access='append') ! solution in concentration
      open(102,file='Arhs.dat',form='formatted',access='append') ! g(y)
      open(103,file='Atmscl.dat',form='formatted',access='append') ! tau
      open(104,file='Asolm.dat',form='formatted',access='append') ! solution in molar concentration
      open(105,file='Aeigen.dat',form='formatted',access='append') ! eigenvalues
      open(106,file='APointers.dat',form='formatted',access='append') ! 
      open(107,file='ATPI.dat',form='formatted',access='append') ! 
      open(108,file='AAPI.dat',form='formatted',access='append') ! 
      open(109,file='AFi.dat',form='formatted',access='append') ! 
      open(110,file='ASS_SyrGlu.dat',form='formatted',access='append')
      open(111,file='AratesLacI.dat',form='formatted',access='append') ! sign of stoich
      open(112,file='ANofExhMod.dat',form='formatted',access='append')
      open(113,file='ADiagnPoin.dat',form='formatted',access='append')
      open(114,file='ASortDiagPP.dat',form='formatted',access='append') ! Sorted diagnostics per point
      open(115,file='AGLCRat.dat',form='formatted',access='append') ! Astrocytes and Neurons Glycolytic Ratios
      open(116,file='AVcapacities.dat',form='formatted',access='append')
      open(117,file='AII.dat',form='formatted',access='append')
      open(118,file='ASlowPoint.dat',form='formatted',access='append')
      open(119,file='ATrend/SS.dat',form='formatted',access='append')
      open(120,file='ATrend/NA.dat',form='formatted',access='append')
      open(1201,file='ATrend/Rat.dat',form='formatted',access='append')    
      open(131,file='AOxidRates.dat',form='formatted',access='append') 
      open(132,file='ARates.dat',form='formatted',access='append')   
c
c     for smart print in different files for every mode
      open(121,file='ADiag/Mod1.dat',form='formatted',access='append')
      open(122,file='ADiag/Mod2.dat',form='formatted',access='append')
      open(123,file='ADiag/Mod3.dat',form='formatted',access='append')
      open(124,file='ADiag/Mod4.dat',form='formatted',access='append')
      open(125,file='ADiag/Mod5.dat',form='formatted',access='append')
      open(126,file='ADiag/Mod6.dat',form='formatted',access='append')
      open(127,file='ADiag/Mod7.dat',form='formatted',access='append')
      open(128,file='ADiag/Mod8.dat',form='formatted',access='append')
      open(129,file='ADiag/Mod9.dat',form='formatted',access='append')
      open(130,file='ADiag/Mod10.dat',form='formatted',access='append')


cpd--------------------------------------------------
cpd   Compute eigenvalues-eigenvectors-timescales
      call jac(n,t,yy,rjac)
      call eigen(rjac,alpha,betta,rlmod,reim)
      do i=1,n
       tau(i,1)=1.0d0/dsqrt(reim(i,1)*reim(i,1)+reim(i,2)*reim(i,2))
       if (reim(i,1).lt.0.0d0) then
        tau(i,2)=-1.0d0
       else
        tau(i,2)=1.0d0
       endif
      enddo
      call fex(n,t,yy,ydot)
      call rates(n,t,yy,k,RR)
c      EPS=qu1/qu2 ! EPSILON
cpd--------------------------------------------------
cpd   Print solution, rates, rhs, eigenvalues and timescales   
      write(101,30) t,yy(1:n)
      write(102,30) t,ydot(1:n)
      write(105,32) t,(reim(i,1:2),i=1,n)
      write(103,30) t,tau(1:n,1)
      write(111,34) t,RR(14)-RR(28),RR(15)-RR(29),RR(16)-RR(30)
      write(131,43) t,RR(8),RR(10),RR(17),RR(18),RR(7),RR(9)  
      write(132,44) t,(RR(i),i=1,k)          
cpd--------------------------------------------------
cpd   Solution in molar concentrations      
      call mol2conc(n,t,yy,yym) 
      write(104,31) t,yym(1:n+2)
cpd--------------------------------------------------                       CSP basis vectors for the criterion
cpd   POinter to feed only pointer needed!!
      do i=1,n
       do j=1,n
        po(i,j)=betta(i,j)*alpha(j,i)  
       enddo
      enddo
cpd   number of exhausted modes with CSP vectors
      CSPrtol=1.0d-2 
      CSPatol=1.0d-16
      call NEM_CSPbv(n,k,t,yy,CSPrtol,CSPatol,rlmod,po,noEM)      

      write(112,35) t,noEM,index 
      po(:,:)=0.0d0
cpd--------------------------------------------------                       CSP basis vectors for the criterion
cpd   Diagnostics from eigenvectors 
      call diag(n,k,NoEM,t,yy,alpha,betta,reim,po,tpi,api,impi,fip)
      temp=1.0d-40
      do i=1,n
       if(temp.lt.dabs(fip(i))) then
        temp=dabs(fip(i))
        index=i
       endif
      enddo
      write(112,35) t,noEM,index      
cpd--------------------------------------------------
cpd   Print diagnostics Pointer, TPI, API, II and Fi 
      call Print_Diag(n,k,t,yy,po,tpi,api,impi,fip,106,107,108,117,109)
cpd--------------------------------------------------
cpd   CSP basis vectors for calculating Slow pointer!!!
      allocate(Ar1(n,NoEM),As1(n,n-NoEM),Br1(NoEM,n),Bs1(n-NoEM,n),
     - Ar2(n,NoEM),As2(n,n-NoEM),Br2(NoEM,n),Bs2(n-NoEM,n)) 
      call CSP_kern(n,k,t,yy,NoEM,PO,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,1) 
cpd   slow Pointer 
      SlowPO(:,:)=0.0d0
      call smult(n,n-NoEM,n,As1,Bs1,SlowPO,n,n-NoEM,n)
      sum1=0.0d0
      do i=1,n
       sum1=sum1+SlowPo(i,i)
      enddo
      SlowPo(:,:)=SlowPo(:,:)/sum1
      deallocate(Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2)   
cpd--------------------------------------------------
cpd Glycolytic Ratios calculation and write
      call GLCRatio(T,n,m,YY,Ras,Rne)
      write(115,41) T,Ras,Rne
cpd--------------------------------------------------
cpd  write the V capacities aVmh,nVmh,aVmL,nVmL
      write(116,*) t,Va,Vn,Val,Vnl
cpd--------------------------------------------------
cpd   START LSODE LOOP
   10 continue 
    
      call DLSODE (FEX, n, yy, t, tout, ITOL, RTOL, ATOL, ITASK,
     1               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
     
      if (ISTATE.lt.0)  go to 98
cpd--------------------------------------------------
cpd   Compute eigenvalues-eigenvectors-timescales
      call jac(n,t,yy,rjac)
      call eigen(rjac,alpha,betta,rlmod,reim)
      do i=1,n
       tau(i,1)=1.0d0/dsqrt(reim(i,1)*reim(i,1)+reim(i,2)*reim(i,2))
       if (reim(i,1).lt.0.0d0) then
        tau(i,2)=-1.0d0
       else
        tau(i,2)=1.0d0
       endif
      enddo
      call fex(n,t,yy,ydot)
      call rates(n,t,yy,k,RR)

cpd--------------------------------------------------
cpd   Print solution, rates, rhs, eigenvalues and timescales   
      write(101,30) t,yy(1:n)
      write(102,30) t,ydot(1:n)
      write(105,32) t,(reim(i,1:2),i=1,n)
      write(103,30) t,tau(1:n,1)
      write(111,34) t,RR(14)-RR(28),RR(15)-RR(29),RR(16)-RR(30)
      write(131,43) t,RR(8),RR(10),RR(17),RR(18),RR(7),RR(9)  
      write(132,44) t,(RR(i),i=1,k)              
cpd
cpd   Solution in molar concentrations      
      call mol2conc(n,t,yy,yym) 
      write(104,31) t,yym(1:n+2)
cpd--------------------------------------------------
cpd   Solution in molar concentrations      
      call mol2conc(n,t,yy,yym) 
      write(104,31) t,yym(1:n+2)
cpd--------------------------------------------------                       CSP basis vectors for the criterion
cpd   POinter to feed only pointer needed!!
      do i=1,n
       do j=1,n
        po(i,j)=betta(i,j)*alpha(j,i)  
       enddo
      enddo
cpd   number of exhausted modes with CSP vectors
      call NEM_CSPbv(n,k,t,yy,CSPrtol,CSPatol,rlmod,po,noEM)    
      write(112,35) t,noEM,index 
cpd--------------------------------------------------                       CSP basis vectors for the criterion      
cpd   Diagnostics from eigenvectors 
      call diag(n,k,NoEM,t,yy,alpha,betta,reim,po,tpi,api,impi,fip)
      temp=1.0d-40
      do i=1,n
       if(temp.lt.dabs(fip(i))) then
        temp=dabs(fip(i))
        index=i
       endif
      enddo
      write(112,35) t,noEM,index        
cpd--------------------------------------------------
cpd   Print diagnostics Pointer, TPI, API, II and Fi 
      call Print_Diag(n,k,t,yy,po,tpi,api,impi,fip,106,107,108,117,109)
cpd--------------------------------------------------
cpd   CSP basis vectors
      allocate(Ar1(n,NoEM),As1(n,n-NoEM),Br1(NoEM,n),Bs1(n-NoEM,n),
     - Ar2(n,NoEM),As2(n,n-NoEM),Br2(NoEM,n),Bs2(n-NoEM,n)) 
      call CSP_kern(n,k,t,yy,NoEM,PO,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,1)
cpd   slow Pointer 
      SlowPO(:,:)=0.0d0
      call smult(n,n-NoEM,n,As1,Bs1,SlowPO,n,n-NoEM,n)
      sum1=0.0d0
      do i=1,n
       sum1=sum1+SlowPo(i,i)
      enddo
      SlowPo(:,:)=SlowPo(:,:)/sum1
      deallocate(Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2)     
cpd--------------------------------------------------
cpd Glycolytic Ratios calculation and write
      call GLCRatio(T,n,m,YY,Ras,Rne)
      write(115,41) T,Ras,Rne 
cpd--------------------------------------------------
cpd  write the V capacities aVmh,nVmh,aVmL,nVmL
      write(116,*) t,Va,Vn,Val,Vnl     
cpd--------------------------------------------------
cpd   Diagnpstics per point-Choose points first
      if (tout.lt.1.0001d+3.and.tout.gt.1.0d+3) go to 200
      if (tout.lt.1030.1d0.and.tout.gt.103.0d+1) go to 200
      if (tout.lt.1050.1d0.and.tout.gt.105.0d+1) go to 200 
      if (tout.lt.1100.1d0.and.tout.gt.11.0d+2) go to 200 
      if (tout.lt.1190.1d0.and.tout.gt.119.0d+1) go to 200          
      if (tout.lt.1280.1d0.and.tout.gt.128.0d+1) go to 200
      if (tout.lt.1330.1d0.and.tout.gt.133.0d+1) go to 200
      if (tout.lt.1450.1d0.and.tout.gt.145.0d+1) go to 200
      if (tout.lt.1500.1d0.and.tout.gt.150.0d+1) go to 200 
      if (tout.lt.1600.1d0.and.tout.gt.160.0d+1) go to 200
      if (tout.lt.1650.1d0.and.tout.gt.165.0d+1) go to 200        
      if (tout.lt.1800.1d0.and.tout.gt.180.0d+1) go to 200                      
      go to 201
cpd   Sort and print the diagnostics in selected time
  200 call sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,PoI,tpiI,
     - apiI,impiI)
      call smart_print(n,k,t,yy,PoS,tpiS,apiS,impiS,PoI,tpiI,
     - apiI,impiI,121,130)

c      call rates(n,t,yy,k,RR)
c      do i=1,k 
c      write(*,*) i, RR(i)
c      enddo    
      write(114,*) 'Point :    ', tout
      write(114,*) 'CSP Pointer Per Species'
      write(114,*) '----------------------------------'
      do i=1,n
       write(114,'(A,i4)') ' Spec :', i
       write(114,37) PoS(i,:)
       write(114,38) PoI(i,:)
      enddo 
      write(114,*) '-------------------------------------------------'
      write(114,*) 'CSP TPI Per Reaction'
      write(114,*) '----------------------------------'
      do i=1,n
       write(114,'(A,i4)') ' Mode :', i
       write(114,39) tpiS(i,:)
       write(114,40) tpiI(i,:)
       write(114,*)
      enddo 
      write(114,*) '-------------------------------------------------'
      write(114,*) 'CSP API Per Reaction'
      write(114,*) '----------------------------------'
      do i=1,n
       write(114,'(A,i4)') ' Mode :', i
       write(114,39) apiS(i,:)
       write(114,40) apiI(i,:)
       write(114,*)
      enddo  
      write(114,*) '-------------------------------------------------'
      write(114,*) 'CSP II Per Reaction'
      write(114,*) '----------------------------------'
      do i=1,n
       write(114,'(A,i4)') ' Spec :', i
       write(114,39) impiS(i,:)
       write(114,40) impiI(i,:)
       write(114,*)
      enddo  
      write(118,*) '-------------------------------------------------'
      write(118,*) '-------------------------------------------------'
      write(118,*) 'Point :    ', tout
      write(118,*) 'CSP Slow Pointer for slow modes'
      write(118,*) '----------------------------------'

      write(118,37) (SlowPO(i,i),i=1,n)


  201 continue 
cpd--------------------------------------------------
cpd   Trends on Lac^e, Lac^a, Lac^n, R^7-R^10 and R^15-R^18 (SS value + in the middle of activation)       
      write(*,*) T
      if (tout.lt.11.0001d+3.and.tout.gt.11.0d+3) go to 300
      if (tout.lt.11110.1d0.and.tout.gt.1111.0d+1) go to 301 
      go to 302   
  300 continue
c      call sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,PoI,tpiI,
c     - apiI,impiI)
      call mol2conc(n,t,yy,yym)
      call rates(n,t,yy,k,RR)
      write(119,42) T, c2, yym(6), yym(8), yym(10), RR(15)-RR(29),
     - RR(16)-RR(30), RR(17), RR(18), impi(10,9), impi(10,18), 
     - impi(10,17), impi(10,11), impi(10,7), impi(10,16), impi(10,30),
     - impi(10,25), api(7,11), api(7,25), api(7,12), api(7,16),
     - api(7,30), api(7,26)
      go to 302
  301 continue
      call sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,PoI,tpiI,
     - apiI,impiI)
      call mol2conc(n,t,yy,yym)
      call rates(n,t,yy,k,RR)
      write(120,42) T, c2, yym(6), yym(8), yym(10), RR(15)-RR(29),
     - RR(16)-RR(30), RR(17), RR(18), impi(10,9), impi(10,18), 
     - impi(10,17), impi(10,11), impi(10,17), impi(10,16), impi(10,30),
     - impi(10,25), api(7,11), api(7,25), api(7,12), api(7,16),
     - api(7,30), api(7,26)
      write(1201,30) T, c2, RR(7), RR(18), RR(17), RR(16), RR(11), 
     - RR(30), RR(9), RR(25), c2
  302 continue


cpd--------------------------------------------------
cpd   Define timestep
      dt=1.0d+1
      if (tout.gt.11.0d+3.and.tout.lt.13.5d+3) DT=1.0d-1
      if (tout.lt.1.0d+2) DT=1.0d-1 ! define time step
      if (tout.lt.1.0d-1) DT=1.0d-3 ! define time step
      if (tout.lt.1.0d-3) DT=1.0d-5 ! define time step
      if (tout.lt.2.0d+3.and.tout.gt.1.0d+3) DT=1.0d-1
cpd--------------------------------------------------
cpd   Define next iteration   
      t=tout
      tout=tout+dt
      if(t.lt.tend) go to 10  ! Define tend above
cpd   END LSODE LOOP
cpd--------------------------------------------------  
      call mol2conc(n,t,yy,yym)
      write(110,*) c1,yym(11),yym(9)     
c
      write(*,*) 'You may now close this window'   
c
cpd------------------------------------------------------------------
cpd   formats for all the .dat
   30 format(E18.11,10('    ',E25.16))   
   31 format(E18.11,12('    ',E25.16))
   32 format(E18.11,20('    ',E25.16))
   33 format(E18.11,a8,i3,30('    ',E25.16))
   34 format(E18.11,3('    ',E25.16))
   35 format(E18.11,2('    ',i3))
   36 format(E18.11,a8,i3,10('    ',E25.16))
   37 format(10(E15.6,'  '))
   38 format(10(i15,'  '))
   39 format(30(E15.6,'  '))
   40 format(30(i15,'  '))
   41 format(E18.11,2('   ',E12.5))
   42 format(E18.11,22('   ',E12.5))   
   43 format(E18.11,6('    ',E25.16)) 
   44 format(E18.11,30('    ',E25.16))     
c   
   99 format(///'t=',D15.8,'    Error halt.. ISTATE =',I3)      
      go to 100
   98 write(*,99)  t, ISTATE
  100 continue   
cpd--------------------------------------------------   
cpd   close all the files
      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      close(116)
      close(117)
      close(118)
      close(119)
      close(120)
      close(121)
      close(122)
      close(123)
      close(124)
      close(125)
      close(126)
      close(127)
      close(128)
      close(129)
      close(130)
      close(131)
      close(132)
    
      STOP
      END
cpd------------------------------------------------------------------
cpd   END of main routine
c
cpd------------------------------------------------------------------
cpd   Subroutine for sorting diagnostics
cpd   Inputs:  n: number of equations
cpd            k: number of reactions
cpd            Po, tpi, api: matrices to be sorted
c
cpd   Outputs: PoS, tpiS, apiS: sorted per rows diagnostic matrices 
cpd            PoI, tpiI, apiI: keeping trace of the reactions. Index matrices  
      subroutine sortDiag(n,k,Po,tpi,api,impi,PoS,tpiS,apiS,impiS,
     - PoI,tpiI,apiI,impiI)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: Po(n,n), tpi(n,k), api(n,k),
     - impi(n,k)
      double precision, intent(out) :: PoS(n,n), tpiS(n,k), apiS(n,k),
     - impiS(n,k)
      integer, intent(out) :: PoI(n,n), tpiI(n,k), apiI(n,k),impiI(n,k)
      integer :: i,j
      double precision :: PoRow(n), tpiRow(k), apiRow(k), impiRow(k)
      integer :: PoInd(n),tpiInd(k),apiInd(k), impiInd(k)
cpd--------------------------------------------------
cpd   Slice the matrices and sort over rows
      do i=1,n      
       PoRow(:)=Po(i,:)
       tpiRow(:)=tpi(i,:)
       apiRow(:)=api(i,:)
       impiRow(:)=impi(i,:)
cpd--------------------------------------------------
cpd    Sort over each slice
       call SortAbs(PoRow,n,PoInd)
       call SortAbs(tpiRow,k,tpiInd)
       call SortAbs(apiRow,k,apiInd)
       call SortAbs(impiRow,k,impiInd)
cpd--------------------------------------------------
cpd    Bind the slices into matrices and their indices        
       PoS(i,:)=PoRow(:)
       PoI(i,:)=PoInd(:)
       tpiS(i,:)=tpiRow(:)
       tpiI(i,:)=tpiInd(:)
       apiS(i,:)=apiRow(:)
       apiI(i,:)=apiInd(:)
       impiS(i,:)=impiRow(:)
       impiI(i,:)=impiInd(:)
      enddo 
cpd-------------------------------------------------- 
      return
      end    
cpd------------------------------------------------------------------      
c
c
cpd------------------------------------------------------------------
cpd   Criterion for the number of exhausted modes.
      subroutine NEM(n,t,yy,CSPrtol,CSPatol,alpha,betta,rlmod,noEM)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: t, yy(n),CSPrtol,CSPatol,
     - alpha(n,n),betta(n,n),rlmod(n)
      integer, intent(out) :: noEM
      integer :: i,j,m1,m2,m
      double precision :: yyer(n),temp(n),fi(n),ydot(n)


      call fex(n,t,yy,ydot)
      call smult(n,n,1,betta,ydot,fi,n,n,1)
      do i=1,n
       yyer(i)=CSPrtol*yy(i)+CSPatol        !Error
      enddo

      noEM=0
  200 continue
      m1=noEM+1
      m2=noEM+2                               !Next two modes
      do i=1,n
       temp(i)=0.0d0
      enddo
      do i=1,m1                               ! evaluate next mode a.f
       do j=1,n
        temp(j)=temp(j)+alpha(j,i)*fi(i)
       enddo
c       write(*,*) alpha(1:n,i) 
c       write(*,*) i, fi(i)
      enddo

      do j=1,n
       temp(j)=temp(j)/dabs(rlmod(m2))              ! multiply with tau+1
      enddo  
c      if(NoEM.eq.8) then
c       do i=1,NoEM
c        write(*,*) temp(i)
c       enddo
c      endif      
c      do i=1,n
c       yyer(i)=dabs(rlmod(m2)/rlmod(m1))*yy(i)+CSPatol        !Error
c      enddo
      M=0
      do j=1,n
       if(dabs(temp(j)).lt.yyer(j)) m=m+1    !|Sum alpha*f|<rtol*y+atol
c       write(*,*) temp(j), yyer(j)
c       write(*,*) j, fi(j)
      enddo

      if(m.eq.n) then
       noEM=noEM+1
       go to 200
      endif

c      write(*,*) NoEM
      return
      end


cpd

cpd   SUBROUT used by lsode to solve the system
cpd   contains the RHS of the ODE system dy/dt=g(y) 
      subroutine FEX(n,t,yy,ydot)
      implicit none
      integer,intent(in) :: n
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: ydot(n)
cpd
      integer :: k,i            
      common k            ! number of reactions passed with common
cpd            
      double precision :: RR(k), st(n,k)

c      include 'paramet.i'      
cpd
cpd   SUB_BLM contains these subs
cpd   generating fex through 
      call rates(n,t,yy,k,RR)
      call stoic(n,k,st)
      call smult(n,k,1,st,RR,ydot,n,k,1)
c      RR(8)=0.0d0
cpd
      return
      end
cpd------------------------------------------------------------------
cpd   SUBROUT used by lsode to solve the system, 
cpd   contains the jacobian of the RHS of the ODE system dy/dt=g(y) 
      subroutine JEX(n,t,yy,ML,MU,pd,NRPD)
      implicit none
      integer, intent(in) :: n, ML, MU, NRPD
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: pd(n,n)
      integer ::i,j
cpd
      integer :: k            
      common k 
      double precision :: st(n,k),gR(k,n),pd1(n,n)
cpd 
cpd      
      do i=1,n
       do j=1,n
        pd(i,j)=0.0D0
       enddo
      enddo
cpd
cpd   SUB_BLM contains these subs
      call stoic(n,k,st)
      call gradR(n,t,yy,k,gR)
      call smult(n,k,n,st,gR,pd,n,k,n)
cpd      
      return
      end	 	  
cpd------------------------------------------------------------------
cpd   SUBROUT contains the jacobian of the RHS of the ODE system dy/dt=g(y) 
cpd   same with jex, user-use because jex is altered after calculations
      subroutine  jac(n,t,yy,pd)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: pd(n,n)
      integer ::i,j
cpd      
      integer :: k            
      common k 
      double precision :: st(n,k),gR(k,n),pd1(n,n)
cpd      
      do i=1,n
       do j=1,n
        pd(i,j)=0.0D0
       enddo
      enddo
cpd
cpd   SUB_BLM contains these subs
      call stoic(n,k,st)
      call gradR(n,t,yy,k,gR)
      call smult(n,k,n,st,gR,pd,n,k,n)
c      write(*,*)'Jac check'
c      write(*,*) pd(1,:)
c      write(*,*) pd(2,:)
c      write(*,*) pd(9,:)
c      write(*,*) pd(10,:)
cpd      
      return
      end
