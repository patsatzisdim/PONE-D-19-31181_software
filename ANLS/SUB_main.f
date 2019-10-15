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
      double precision, dimension(n,n) :: pd,rjac,alpha,betta,test,po
      integer :: ML, MU, NRPD
      double precision, dimension(n,2) :: reim,tau
      double precision :: tpi(n,m),api(n,m),fip(n),RR(m)
      character :: varName*3
      double precision :: CSPrtol,CSPatol,temp
      integer :: noEM,index
      double precision :: PoS(n,n), tpiS(n,m), apiS(n,m)
      integer :: PoI(n,n), tpiI(n,m), apiI(n,m)
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
cpd
cpd   Solution in molar concentrations      
      call mol2conc(n,t,yy,yym) 
      write(104,31) t,yym(1:n+2)
cpd
cpd   CSP Diagnostics
      call diag(n,k,t,yy,alpha,betta,reim,po,tpi,api,fip) 
cpd--------------------------------------------------
cpd   Print diagnostics Pointer, TPI, API and Fi 
c      write(106,30) t, (po(i,i),i=1,n)              ! Diagonal Pointers
      do i=1,n
       write(106,36) t, 'Spec  ',i, po(i,1:n) 
       write(107,33) t, 'Spec  ',i, tpi(i,1:k)
       write(108,33) t, 'Spec  ',i, api(i,1:k)
      enddo
      write(109,30) t, fip(1:n)
cpd
cpd   CSP Exhausted modes criterion
      CSPrtol=1.0d-1 
      CSPatol=1.0d-20
      call NEM(n,t,yy,CSPrtol,CSPatol,alpha,fip,rlmod,noEM) 
      temp=1.0d-40
      do i=1,n
       if(temp.lt.fip(i)) then
        temp=fip(i)
        index=i
       endif
      enddo
      write(112,35) t,noEM,index     
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
cpd   Change the parameters.
c
c
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
cpd
cpd   Solution in molar concentrations      
      call mol2conc(n,t,yy,yym) 
      write(104,31) t,yym(1:n+2)
cpd
cpd   CSP Diagnostics
      call diag(n,k,t,yy,alpha,betta,reim,po,tpi,api,fip) 
cpd--------------------------------------------------
cpd   Print diagnostics Pointer, TPI, API and Fi 
c      write(106,30) t, (po(i,i),i=1,n)              ! Diagonal Pointers
      do i=1,n
       write(106,36) t, 'Spec  ',i, po(i,1:n) 
       write(107,33) t, 'Spec  ',i, tpi(i,1:k)
       write(108,33) t, 'Spec  ',i, api(i,1:k)
      enddo
      write(109,30) t, fip(1:n)
cpd
cpd   CSP Exhausted modes criterion
      call NEM(n,t,yy,CSPrtol,CSPatol,alpha,fip,rlmod,noEM) 
      temp=1.0d-40
      do i=1,n
       if(temp.lt.fip(i)) then
        temp=fip(i)
        index=i
       endif
      enddo
      write(112,35) t,noEM,index
cpd--------------------------------------------------
cpd Glycolytic Ratios calculation and write
      call GLCRatio(T,n,m,YY,Ras,Rne)
      write(115,41) T,Ras,Rne 
cpd--------------------------------------------------
cpd  write the V capacities aVmh,nVmh,aVmL,nVmL
      write(116,*) t,Va,Vn,Val,Vnl     
cpd--------------------------------------------------
cpd   Diagnpstics per point-Choose points first
      if (tout.lt.9.51d+2.and.tout.gt.9.50d+2) go to 200
      if (tout.lt.1012.6d0.and.tout.gt.1012.5d0) go to 200
      if (tout.lt.1050.1d0.and.tout.gt.1050.0d0) go to 200
      if (tout.lt.1100.1d0.and.tout.gt.1100.0d0) go to 200
      go to 201
cpd   Sort and print the diagnostics in selected time
  200 call sortDiag(n,k,Po,tpi,api,PoS,tpiS,apiS,PoI,tpiI,apiI)      
      write(114,*) 'Point :    ', tout
      write(114,*) 'CSP Pointer Per Species'
      write(114,*) '----------------------------------'
      do i=1,n
       write(114,'(A,i4)') ' Mode :', i
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
      write(114,*) '-------------------------------------------------'
  201 continue      
cpd--------------------------------------------------
cpd   Define timestep
      dt=1.0d0
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
cpd    perturbation on parameters continue  
c      newVal=0.5d-3     
c      call ChangeParamet(95,'GCF',newVal)
c      
c      if(t.le.tend2) go to 10

      call mol2conc(n,t,yy,yym)
      write(110,*) c1,yym(11),yym(9)     


      write(*,*) 'You may now close this window'   


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
      subroutine sortDiag(n,k,Po,tpi,api,PoS,tpiS,apiS,PoI,tpiI,apiI)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: Po(n,n), tpi(n,k), api(n,k)
      double precision, intent(out) :: PoS(n,n), tpiS(n,k), apiS(n,k)
      integer, intent(out) :: PoI(n,n), tpiI(n,k), apiI(n,k)
      integer :: i,j
      double precision :: PoRow(n), tpiRow(k), apiRow(k), PoRowS(n)
      integer :: PoInd(n),tpiInd(k),apiInd(k)

cpd   Slice the matrices and sort over rows
      do i=1,n      
       do j=1,n
        PoRow(j)=Po(i,j)
       enddo
       do j=1,k
        tpiRow(j)=tpi(i,j)
        apiRow(j)=api(i,j)
       enddo
c
c
cpd    Sort over each slice
       call SortAbs(PoRow,n,PoInd)
       call SortAbs(tpiRow,k,tpiInd)
       call SortAbs(apiRow,k,apiInd)
c
cpd    Bind the slices into matrices and their indices        
       do j=1,n
        PoS(i,j)=PoRow(j)
        PoI(i,j)=PoInd(j)
       enddo
       do j=1,k
        tpiS(i,j)=tpiRow(j)
        tpiI(i,j)=tpiInd(j)
        apiS(i,j)=apiRow(j)
        apiI(i,j)=apiInd(j)
       enddo
      enddo 
c  
      return
      end    
c
c
cpd------------------------------------------------------------------
cpd   Criterion for the number of exhausted modes.
      subroutine NEM(n,t,yy,CSPrtol,CSPatol,alpha,fi,rlmod,noEM)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: t, yy(n),CSPrtol,CSPatol,
     - alpha(n,n),fi(n),rlmod(n)
      integer, intent(out) :: noEM
      integer :: i,j,m1,m2,m
      double precision :: yyer(n),temp(n)


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
      do i=1,m1                               ! evaluate next mode a f
       do j=1,n
        temp(j)=temp(j)+alpha(j,i)*fi(j)
       enddo
      enddo
      do j=1,n
       temp(j)=temp(j)/rlmod(m2)              ! multiply with tau+1
      enddo  
      M=0
      do j=1,n
       if(dabs(temp(j)).lt.yyer(j)) m=m+1    !|Sum alpha*f|<rtol*y+atol
      enddo
      if(m.eq.n) then
       noEM=noEM+1
       go to 200
      endif

      return
      end

cpd------------------------------------------------------------------
cpd   Diagnostics subroutine. Printings are passed out at the main code.
cpd    ----------------------------------------------------------------CODE CHECKED. CASE OF complex EIG is not considered
      subroutine Diag(n,k,t,yy,alpha,betta,reim,po,tpi,bsr,fip)
      implicit none
      integer, intent(in) :: n,k
      integer :: i,j,ii,count,countSP,count1
      double precision, intent(in) :: t, yy(n), reim(n,2)
      double precision, intent(inout) ::  alpha(n,n), betta(n,n)
      double precision, intent(out) :: po(n,n),tpi(n,k),bsr(n,k),fip(n)
      double precision :: dpn(n), st(n,k), gR(k,n), gradSRk(n,n),
     - tempS(n,1), tempGR(1,n), tempB(1,n), BgSRk(1,n), sum1, 
     - BgSRkA(n,k), tsum1, tsum2, sum2(n), RR(k), FI(4,n),
     - tempSR(n,k), dumBSR(n,k)
      
      count=0
      count1=0

cpd--------------------------------------------------
cpd   CSP Pointer through eigenvectors
cpd   Watch out for complex eigenvalues
      do i=1,n
       if (reim(i,2).eq.0.0d0) then     ! REAL EIGENVALUES
        do j=1,n
         po(i,j)=betta(i,j)*alpha(j,i)  ! Pointer through eigenvectors
         if(i.eq.j) dpn(i)=po(i,j)      ! Diagonal element
        enddo
       else                             ! COMPLEX EIGENVALUES
       	count=count+1                   ! counter to check if this is the first complex eigen
        if(mod(count,2).ne.0) go to 60  ! if this is the first complex eigen from the pair, just pass, if not continue
        do j=1,n
         po(i-1,j)=0.5d0*(betta(i-1,j)*alpha(j,i-1)+
     -                    betta(i,j)*alpha(j,i))
         po(i,j)=0.5d0*(betta(i,j)*alpha(j,i-1)-
     -                  betta(i-1,j)*alpha(j,i))
        enddo
   60   continue
       endif                             
      enddo
cpd--------------------------------------------------
cpd
cpd   CSP TPI through eigencvectors
      call stoic(n,k,st)
      call gradR(n,t,yy,k,gR)
cpd
cpd   Loop over the reactions
      do j=1,k
cpd   Construct the grad(s_i R^i)
cpd   First s_i 
       do i=1,n
        tempS(i,1)=st(i,j)
       enddo       
cpd   Second grad(R^i)
       do i=1,n
        tempGR(1,i)=gR(j,i)
       enddo
       call smult(n,1,n,tempS,tempGR,gradSRk,n,1,n)  !  grad(S_i R^i)
      countSP=0
cpd 
cpd   Loop over the species          
       do i=1,n
       	countSP=countSP+1                             ! counter over species        	
       	if (reim(i,2).eq.0.0d0) then                  ! REAL EIGENVALUES
         do ii=1,n
          tempB(1,ii)=betta(countSP,ii)
         enddo
         call smult(1,n,n,tempB,gradSRk,BgSRk,1,n,n)  ! betta^j S_i gradR^i
         sum1=0.0d0  
         do ii=1,n
          sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP)
         enddo
         BgSRkA(countSP,j)=sum1                       ! betta^j S_i gradR^i alpha_j= c^n_k
        else                                          ! COMPLEX EIGENVALUES
         count=count+1                                ! counter to check if this is the first complex eigen
         tsum1=0.0d0
         tsum2=0.0d0                 
         if(mod(count,2).ne.0) go to 70               ! if this is the first complex eigen from the pair, just pass, if not continue         
          do ii=1,n
           tempB(1,ii)=betta(countSP-1,ii)
          enddo
          call smult(1,n,n,tempB,gradSRk,BgSRk,1,n,n)
          sum1=0.0D0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP-1) !b^1 J a_1
          enddo
          tsum1=tsum1+sum1
          sum1=0.0d0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP)   !b^1 J a_2
          enddo
          tsum2=tsum2+sum1
          
          do ii=1,n
           tempB(1,ii)=betta(countSP,ii)
          enddo
          call smult(1,n,n,tempB,gradSRk,BgSRk,1,n,n)
          sum1=0.0D0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP-1) !b^2 J a_1
          enddo
          tsum2=tsum2-sum1							 !b^1 J a_2 - b^2 J a_1
          sum1=0.0D0
          do ii=1,n
           sum1=sum1+BgSRk(1,ii)*alpha(ii,countSP)   !b^2 J a_2
          enddo
          tsum1=tsum1+sum1                           !b^1 J a_1 + b^2 J a_2
          BgSRkA(countSP,j)=0.5d0*tsum1              ! betta^j S_i gradR^i alpha_j= c^1_k
          BgSRkA(countSP+1,j)=0.5d0*tsum2            ! betta^j S_i gradR^i alpha_j= c^2_k
   70    continue
        endif	
        if (countSP.eq.n) go to 80
       enddo
cpd   End of loop over the species
   80 continue
      enddo
cpd   End of loop over the reactions   	
cpd   Now continue in construction of TPI
      
      do i=1,n
       sum2(i)=0.0d0
       do j=1,k
        sum2(i)=sum2(i)+dabs(BgSRkA(i,j))
       enddo
      enddo
      do i=1,n
       do j=1,k
       	tpi(i,j)=BgSRkA(i,j)/sum2(i)
       enddo
      enddo
cpd--------------------------------------------------
cpd            
cpd   API construction through eigenvectors 
cpd   Amplitude construction
      call stoic(n,k,st)
      call rates(n,t,yy,k,RR)

      do i=1,n
       do j=1,k
       	tempSR(i,j)=st(i,j)*RR(j)     ! matrix of s^j_i R^j
       enddo
      enddo
      call smult(n,n,k,betta,tempSR,BSR,n,n,k)    ! b^i s^j_i R^j  
cpd   Amplitudes to be positive
      do i=1,n                       ! Loop over species
       FI(1,i)=0.0d0
       FI(2,i)=0.0d0
       FI(3,i)=0.0d0
       FI(4,i)=0.0d0
       do j=1,k                            
       	FI(1,i)=FI(1,i)+BSR(i,j)             ! Sum over the reactions
       	FI(2,i)=FI(2,i)+dabs(BSR(i,j))       ! abs
       	if(BSR(i,j).gt.0.0d0) FI(3,i)=FI(3,i)+BSR(i,j)      ! positive
       	if(BSR(i,j).lt.0.0d0) FI(4,i)=FI(4,i)+BSR(i,j)      ! negative
       enddo
       sum1=1.0d0/FI(2,i)                       ! 1/sum(abs)
       do j=1,k
       	dumBSR(i,j)=BSR(i,j)                       ! dummy to save the values
        BSR(i,j)=BSR(i,j)*sum1
       enddo	
      enddo   						 ! End of loop over species
      do i=1,n                       ! New loop over species after calculations
       do j=1,k
        BSR(i,j)=BSR(i,j)*dsign(1.d0,FI(1,i))
        dumBSR(i,j)=dumBSR(i,j)*dsign(1.d0,FI(1,i))   ! take the sign of the sum of BSR
       enddo
       do j=1,n
       	alpha(i,j)=alpha(i,j)*dsign(1.d0,FI(1,i))     ! change the sign of eigenvectors
        betta(i,j)=betta(i,j)*dsign(1.d0,FI(1,i))
       enddo
       if(FI(1,i).lt.0.0d0) then
       	sum1=FI(3,i)
		    FI(3,i)=-FI(4,i)
		    FI(4,i)=-sum1
       endif
	     FI(1,i)=FI(1,i)*dsign(1.d0,FI(1,i))            ! F^i positive
       fip(i)=FI(1,i)
      enddo                          ! End of loop over species
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
