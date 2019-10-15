      subroutine CSP_kern(n,k,t,y,Ms,PO,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,
     - NoPh)
      implicit none
      integer, intent(in) :: n,k,Ms,NoPh
      double precision, intent(in) :: t, y(n), PO(n,n)
      double precision :: Ar0(n,Ms),As0(n,n-Ms),Br0(Ms,n),Bs0(n-Ms,n),
     - dBrdt(Ms,n)
      double precision :: Ar(n,Ms),As(n,n-Ms),Br(Ms,n),Bs(n-Ms,n)
      double precision, intent(out) :: Ar1(n,Ms),As1(n,n-Ms),Br1(Ms,n),
     - Bs1(n-Ms,n)
      double precision, intent(out) :: Ar2(n,Ms),As2(n,n-Ms),Br2(Ms,n),
     - Bs2(n-Ms,n)
      integer :: i
      logical :: Br1call

cpd-----------------------------------------------  
c      Prelimineries of CSP: initial Basis Vectors
cpd-----------------------------------------------  
      call Init_BV(n,Ms,PO,Ar0,As0,Br0,Bs0)


c      write(*,*) 'Initial BV',Br0(5,:)
cpd-----------------------------------------------  
c     Phase 1: dBr0/dt=0     Br and Ar refs
cpd-----------------------------------------------  
c     Step 1: B^r refinement
      Br1call=.true.
      call Br_ref(n,k,Ms,t,y,Br1call,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,
     - dBrdt)    
c     Step 2: A^r refinement
      call Ar_ref(n,k,Ms,t,y,Ar,As,Br,Bs,Ar1,As1,Br1,Bs1)
      if(NoPh.eq.1) go to 10
cpd-----------------------------------------------  
c     Phase 2: dBr0/dt=/0    Br ref ONLY!!!
cpd-----------------------------------------------  
c     Step 1: B^r refinement
      Br1call=.false.
      call Br_ref(n,k,Ms,t,y,Br1call,Ar1,As1,Br1,Bs1,Ar2,As2,Br2,Bs2,
     - dBrdt)
c      write(*,*) matmul(Br2,Ar2)
c      write(*,*)
c      write(*,*) matmul(Br2,As2)
c      write(*,*)
c      write(*,*) matmul(Bs2,Ar2)
c      write(*,*)
c      write(*,*) matmul(Bs2,As2)
c      write(*,*)
c      write(*,*) matmul(Ar2,Br2)+matmul(As2,Bs2)
c      stop
c      stop
cpd-----------------------------------------------  
c     Cut the small values in Basis Vectors   !LEAVE THIS FOR NOW
cpd----------------------------------------------- 
   10 continue     
      return 
      end

cpd-------------------------------------------------------------------  
c
c     Subrout for Br refinement 
c
cpd-------------------------------------------------------------------
      subroutine Br_ref(n,k,m,t,y,Br1call,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,
     - dBrdt)
      implicit none	
      integer, intent(in) :: n,m,k
      double precision, intent(in) :: t,y(n),Ar0(n,m),As0(n,n-m),
     - Br0(m,n),Bs0(n-m,n)
      double precision, intent(out) :: Ar(n,m),As(n,n-m),Br(m,n),
     - Bs(n-m,n)
      double precision, intent(inout) :: dBrdt(m,n)   !Output at phase 1, Input at phase 2
      logical, intent(in) :: Br1call
      double precision :: st(n,k),gR(k,n),djac(n,n),dBr(m,n),Lrr(m,m),
     - Trr(m,m),dJacdt(n,n)
      double precision :: temp1(m,n),Un(n,n),temp2(n,n),temp3(m,n)
      integer :: IPVT(n),INFO
cpd-----------------------------------------------
cpd   Calculate Lamba^r_r
      call stoic(n,k,st)
      call gradR(n,t,y,k,gR)
      call smult(n,k,n,st,gR,djac,n,k,n)      ! Jacobian           
cpd-----------------------------------------------
      if (Br1call) then
       dBr(:,:)=0.0d0
      else
       dBr(:,:)=dBrdt(:,:)
      endif                               ! dBr0/dt
cpd-----------------------------------------------
      call smult(m,n,n,Br0,djac,temp1,m,n,n)  ! Br0.J
      temp1(:,:)=temp1(:,:)+dBr(:,:)          ! Br0.J+dBr0/dt
      call smult(m,n,m,temp1,Ar0,Lrr,m,n,m)   ! (Br0.J+dBr0/dt).Ar0=Lrr    
      call sinve(m,Lrr,Trr,m,IPVT,INFO)       ! Trr=(Lrr)^-1 
c      if (T.gt.1348.9d0.and.T.lt.1349.0d0) write(*,*) Lrr           
cpd----------------------------------------------- 
      call smult(m,m,n,Trr,temp1,Br,m,m,n)    ! Br=Trr.Br0.J
c      call smult(m,n,n,temp4,djac,Br,m,n,n)    
      Ar(:,:)=Ar0(:,:)					      ! Ar=Ar0
      Bs(:,:)=Bs0(:,:)					      ! Bs=Bs0   
      call unitary(n,Un,n)                    
      call smult(n,m,n,Ar,Br,temp2,n,m,n)	  ! Ar.Br
      temp2(:,:)=Un(:,:)-temp2(:,:)           ! I-Ar.Br
      call smult(n,n,n-m,temp2,As0,As,n,n,n-m)! (I-Ar.Br).As0 
cpd-----------------------------------------------
      if(Br1call) then
       temp1(:,:)=0.0d0
       call smult(m,m,n,Trr,Br0,temp1,m,m,n)  !Trr.Br0
       call Djac_dt(n,t,y,dJacdt)
       call smult(m,n,n,temp1,dJacdt,temp3,m,n,n) ! (Trr.Br0).dJ/dt
       call smult(m,n,n,temp3,temp2,dBrdt,m,n,n)       ! (Trr.Br0).dJ/dt.(I-Ar.Br)
      endif
cpd-----------------------------------------------
      return
      end
c
c
cpd-------------------------------------------------------------------  
c
c     Subrout for Ar refinement 
c
cpd-------------------------------------------------------------------
      subroutine Ar_ref(n,k,m,t,y,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs)
      implicit none	
      integer, intent(in) :: n,m,k
      double precision, intent(in) :: t,y(n),Ar0(n,m),As0(n,n-m),
     - Br0(m,n),Bs0(n-m,n)
      double precision, intent(out) :: Ar(n,m),As(n,n-m),Br(m,n),
     - Bs(n-m,n)
      double precision :: st(n,k),gR(k,n),djac(n,n),dBr(m,n),Lrr(m,m),
     - Trr(m,m),dJacdt(n,n)
      double precision :: temp1(m,n),Un(n,n),temp2(n,m),temp3(n,n)
      integer :: IPVT(n),INFO
cpd-----------------------------------------------
cpd   Calculate Lamba^r_r
      call stoic(n,k,st)
      call gradR(n,t,y,k,gR)
      call smult(n,k,n,st,gR,djac,n,k,n)      ! Jacobian              
cpd-----------------------------------------------
      call smult(m,n,n,Br0,djac,temp1,m,n,n)  ! Br0.J
      call smult(m,n,m,temp1,Ar0,Lrr,m,n,m)   ! (Br0.J).Ar0=Lrr
      call sinve(m,Lrr,Trr,m,IPVT,INFO)       ! Trr=(Lrr)^-1     
cpd----------------------------------------------- 

      call smult(n,n,m,djac,Ar0,temp2,n,n,m)   
      call smult(n,m,m,temp2,Trr,Ar,n,m,m)     ! Ar=(J.Ar0).Trr
      As(:,:)=As0(:,:)					       ! As=As0
      Br(:,:)=Br0(:,:)					       ! Br=Br0   
      call unitary(n,Un,n)                    
      call smult(n,m,n,Ar,Br,temp3,n,m,n)	   ! Ar.Br
      temp3(:,:)=Un(:,:)-temp3(:,:)            ! I-Ar.Br
      call smult(n-m,n,n,Bs0,temp3,Bs,n-m,n,n) ! Bs0.(I-Ar.Br) 

c      call smult(m,m,n,Trr,temp1,Br,m,m,n)    ! Br=Trr.(Br0.J+dBr0/dt)
c      Ar(:,:)=Ar0(:,:)					       ! Ar=Ar0
c      Bs(:,:)=Bs0(:,:)					       ! Bs=Bs0   
c      call unitary(n,Un,n)                    
c      call smult(n,m,n,Ar,Br,temp2,n,m,n)	   ! Ar.Br
c      temp2(:,:)=Un(:,:)-temp2(:,:)           ! I-Ar.Br
c      call smult(n,n,n-m,temp2,As0,As,n,n,n-m)! (I-Ar.Br).As0 
cpd-----------------------------------------------
      return
      end
c
c      
cpd-------------------------------------------------------------------  
c
c     Subrout for initial set of Basis vectors. 
c     According to the greatest specie of the PO in each mode.
c
cpd------------------------------------------------------------------- 
      subroutine Init_BV(n,m,PO,Ar,As,Br,Bs)
      implicit none
      integer, intent(in) :: n,m
      double precision, intent(in) :: PO(n,n)
      double precision, intent(out) :: Ar(n,m),As(n,n-m),Br(m,n),
     - Bs(n-m,n)
      integer :: i,maxLocPo
cpd-----------------------------------------------  

      Ar(1:n,1:m)=0.0d0
      As(1:n,1:n-m)=0.0d0
      Br(1:m,1:n)=0.0d0
      Bs(1:n-m,1:n)=0.0d0
cpd-----------------------------------------------        
cpd   read CSP pointer and put 1's and 0's for the species of every mode
      do i=1,m
       maxLocPo=maxloc(PO(i,:),n)          ! the greater           
       Ar(maxLocPo,i)=1.0d0                 
       Br(i,maxLocPo)=1.0d0
      enddo
      do i=1,n-m
       maxLocPo=maxloc(PO(i+m,:),n)          ! the greater           
       As(maxLocPo,i)=1.0d0                 
       Bs(i,maxLocPo)=1.0d0
      enddo
cpd-----------------------------------------------  
      return
      end
cpd-----------------------------------------------     
c
c      
cpd-------------------------------------------------------------------  
c
c     Subrout for the calculation of the number of exhausted 
c     modes through the CSP basis vectors 
c
cpd------------------------------------------------------------------- 
      subroutine NEM_CSPbv(n,k,t,yy,CSPrtol,CSPatol,rlmod,po,noEM)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t, yy(n),CSPrtol,CSPatol,
     - po(n,n),rlmod(n)
      integer, intent(out) :: noEM
      integer :: i,j,m1,m2,m,NoPh
      double precision :: yyer(n),temp(n),fi(n),ydot(n)
      double precision, dimension(:,:), allocatable :: Ar0,As0,Br0,Bs0,
     -Ar,As,Br,Bs,dummy,Ar1,As1,Br1,Bs1
      double precision, dimension(:), allocatable :: fr

c      do i=1,n
c       yyer(i)=CSPrtol*yy(i)+CSPatol        !Error
c      enddo
cpd-----------------------------------------------


c      call smult(n,n,1,betta,ydot,fi,n,n,1)

cpd-----------------------------------------------
      noEM=0
  200 continue
      m1=noEM+1
      m2=noEM+2                               !Next two modes
      do i=1,n
       temp(i)=0.0d0
      enddo
cpd-----------------------------------------------      
c     calculate CSP basis vectors (1-ref) with m1 exhausted modes
      allocate(Ar0(n,m1),As0(n,n-m1),Br0(m1,n),Bs0(n-m1,n),Ar(n,m1),
     -As(n,n-m1),Br(m1,n),Bs(n-m1,n),dummy(m1,n),Ar1(n,m1),As1(n,n-m1),
     -Br1(m1,n),Bs1(n-m1,n),fr(m1))
c      
      call Init_BV(n,m1,po,Ar0,As0,Br0,Bs0)
c      Br1call=.true.
      dummy(:,:)=0.0d0
c      write(*,*) T
      call Br_ref(n,k,m1,t,yy,.false.,Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,
     - dummy)
      call Ar_ref(n,k,m1,t,yy,Ar,As,Br,Bs,Ar1,As1,Br1,Bs1)          ! CSP bv
      call fex(n,t,yy,ydot)
      call smult(m1,N,1,Br1,ydot,fr,m1,N,1)                         ! fr
      call smult(n,m1,1,Ar1,fr,temp,n,m1,1)                         ! ar.fr
      deallocate(Ar0,As0,Br0,Bs0,Ar,As,Br,Bs,dummy,Ar1,As1,Br1,Bs1,fr) 

cpd-----------------------------------------------

      do j=1,n
       temp(j)=temp(j)/dabs(rlmod(m2))              ! multiply with tau+1
      enddo  
c      if(NoEM.eq.8) then
c       do i=1,NoEM
c        write(*,*) temp(i)
c       enddo
c      endif      
      do i=1,n
       yyer(i)=3.0d-1*dabs(rlmod(m2)/rlmod(m1))*yy(i)+CSPatol        !Error
      enddo
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



