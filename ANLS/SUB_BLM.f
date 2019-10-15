cpd-----------This subs contain the problem--------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used to generate the Rates
      subroutine rates(n,t,yy,k,RR)
      implicit none 
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: RR(k)
      double precision :: y1,y2,y3,y4,y5,y6,y7,y8,y9,y10
      
      include 'paramet.i'

cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
      y5=yy(5)
      y6=yy(6)
      y7=yy(7)
      y8=yy(8) 
      y9=yy(9)
      y10=yy(10)  
cpd
cpd   Reaction Rates from mathematica
      RR(1)=(c1*(k1*V1 + y1))/(c1*cr2*k1*V1 + cr1*k1**2*V1 + 
     - c1*cr4*y1 + cr3*k1*y1)  
      RR(2)=(y1*(k1*V2 + y2))/(cr1*k1**2*V1*V2 + cr3*k1*V2*y1 +
     - cr2*k1*V1*y2 + cr4*y1*y2)
      RR(3)=(y2*(k2*V3 + y3))/(cr5*k2**2*V2*V3 + cr6*k2*V3*y2 +
     - cr7*k2*V2*y3 + cr8*y2*y3)
      RR(4)=(k3*y2)/(V2*1.0d+15)
      RR(5)=(y3*(k2*V4 + y4))/(cr9*k2**2*V3*V4 + cr11*k2*V4*y3 +
     - cr10*k2*V3*y4 + cr12*y3*y4)
      RR(6)=(y4*(k4*V5 + y5))/(cr13*k4**2*V4*V5 + cr14*k4*V5*y4 +
     - cr15*k4*V4*y5 + cr16*y4*y5)
      RR(7)=(Va*y3)/(k5*V3 + y3)
      RR(8)=fac*(Va*y3)/(12.0d0*(k5*V3 + y3))
      RR(9)=(Vn*y5)/(k5*V5 + y5)
      RR(10)=(Vn*y5)/(12.0d0*(k5*V5 + y5))
      RR(11)=(c2*(k6*V1 + y6))/(k6*(c2*cr18 + cr17*k6)*V1 +
     - (c2*cr20 + cr19*k6)*y6)
      RR(12)=(y6*(k6*V2 + y7))/(k6*V2*(cr17*k6*V1 + cr19*y6) +
     - (cr18*k6*V1 + cr20*y6)*y7)
      RR(13)=(y7*(k7*V3 + y8))/(k7*V3*(cr21*k7*V2 + cr22*y7) +
     - (cr23*k7*V2 + cr24*y7)*y8)
      RR(14)=(k3*y7)/(V2*1.0d+15)
      RR(15)=(y8*(k7*V4 + y9))/(k7*V4*(cr25*k7*V3 + cr27*y8) +
     - (cr26*k7*V3 + cr28*y8)*y9)
      RR(16)=((k8*V5 + y10)*y9)/(k8*V4*(cr29*k8*V5 +
     - cr31*y10) + (cr30*k8*V5 + cr32*y10)*y9)
      RR(17)=(Val*y8)/(k9*V3 + y8)
      RR(18)=(Vnl*y10)/(k9*V5 + y10)      
      RR(19)=((c1 + k1)*y1)/(c1*cr2*k1*V1 + cr1*k1**2*V1 + c1*cr4*y1 +
     - cr3*k1*y1)
      RR(20)=((k1*V1 + y1)*y2)/(cr1*k1**2*V1*V2 + cr3*k1*V2*y1 +
     - cr2*k1*V1*y2 + cr4*y1*y2)
      RR(21)=((k2*V2 + y2)*y3)/(cr5*k2**2*V2*V3 + cr6*k2*V3*y2 +
     - cr7*k2*V2*y3 + cr8*y2*y3)
      RR(22)=(k3*y4)/(V4*1.0d+15)
      RR(23)=((k2*V3 + y3)*y4)/(cr9*k2**2*V3*V4 + cr11*k2*V4*y3 +
     - cr10*k2*V3*y4 + cr12*y3*y4)
      RR(24)=((k4*V4 + y4)*y5)/(cr13*k4**2*V4*V5 + cr14*k4*V5*y4 +
     - cr15*k4*V4*y5 + cr16*y4*y5)
      RR(25)=((c2 + k6)*y6)/(k6*(c2*cr18 + cr17*k6)*V1 +
     - (c2*cr20 + cr19*k6)*y6)
      RR(26)=((k6*V1 + y6)*y7)/(k6*V2*(cr17*k6*V1 + cr19*y6) +
     - (cr18*k6*V1 + cr20*y6)*y7)
      RR(27)=((k7*V2 + y7)*y8)/(k7*V3*(cr21*k7*V2 + cr22*y7) +
     - (cr23*k7*V2 + cr24*y7)*y8)
      RR(28)=(k3*y9)/(V4*1.0d+15)
      RR(29)=((k7*V3 + y8)*y9)/(k7*V4*(cr25*k7*V3 + cr27*y8) +
     - (cr26*k7*V3 + cr28*y8)*y9)
      RR(30)=(y10*(k8*V4 + y9))/(k8*V4*(cr29*k8*V5 +
     - cr31*y10) + (cr30*k8*V5 + cr32*y10)*y9)
cpd          
      return
      end
cpd
cpd------------------------------------------------------------------
cpd   SUBROUT used for stoichiometry
      subroutine stoic(n,k,st)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(out) :: st(n,k)
      integer :: i,j
cpd
cpd   initialize first      
      do i=1,n
       do j=1,k
        st(i,j)=0.0d0 
       enddo
      enddo
cpd
cpd   values from mathematica
      st(1,1)=+1.0d0
      st(1,2)=-1.0d0
      st(1,19)=-1.0d0
      st(1,20)=+1.0d0      
      st(2,2)=+1.0d0
      st(2,3)=-1.0d0
      st(2,4)=-1.0d0
      st(2,20)=-1.0d0
      st(2,21)=+1.0d0
      st(2,22)=+1.0d0
      st(3,3)=+1.0d0
      st(3,5)=-1.0d0
      st(3,7)=-1.0d0
      st(3,8)=-1.0d0
      st(3,21)=-1.0d0 
      st(3,23)=+1.0d0
      st(4,4)=+1.0d0
      st(4,5)=+1.0d0
      st(4,6)=-1.0d0
      st(4,22)=-1.0d0
      st(4,23)=-1.0d0
      st(4,24)=+1.0d0
      st(5,6)=+1.0d0
      st(5,9)=-1.0d0
      st(5,10)=-1.0d0
      st(5,24)=-1.0d0
      st(6,11)=+1.0d0
      st(6,12)=-1.0d0
      st(6,25)=-1.0d0
      st(6,26)=+1.0d0
      st(7,12)=+1.0d0
      st(7,13)=-1.0d0
      st(7,14)=-1.0d0
      st(7,26)=-1.0d0
      st(7,27)=+1.0d0
      st(7,28)=+1.0d0
      st(8,7)=+2.0d0
      st(8,13)=+1.0d0
      st(8,15)=-1.0d0
      st(8,17)=-1.0d0
      st(8,27)=-1.0d0
      st(8,29)=+1.0d0
      st(9,14)=+1.0d0
      st(9,15)=+1.0d0
      st(9,16)=-1.0d0
      st(9,28)=-1.0d0
      st(9,29)=-1.0d0
      st(9,30)=+1.0d0
      st(10,9)=+2.0d0
      st(10,16)=+1.0d0
      st(10,18)=-1.0d0
      st(10,30)=-1.0d0
    
cpd      
      return
      end      
cpd------------------------------------------------------------------
cpd
cpd------------------------------------------------------------------
cpd   SUBROUT used to calulate grad(R) analytically from mathematica
      subroutine gradR(n,t,yy,k,gR)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: gR(k,n)
      double precision :: y1,y2,y3,y4,y5,y6,y7,y8,y9,y10
      integer :: i,j
      
      include 'paramet.i'

cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
      y5=yy(5)
      y6=yy(6)
      y7=yy(7)
      y8=yy(8) 
      y9=yy(9)
      y10=yy(10) 
cpd
      do i=1,k
       do j=1,n
        gR(i,j)=0.0D0
       enddo
      enddo
cpd   Reaction Rates from mathematica     
cpd   gR(#Reaction,#Species)
cpd
      gR(1,1)=(c1*k1*(c1*(cr2 - cr4) + (cr1 - cr3)*k1)*V1)/
     -   (c1*cr2*k1*V1 + cr1*k1**2*V1 + c1*cr4*y1 + cr3*k1*y1)**2
      gR(2,1)=(k1*V1*(k1*V2 + y2)*(cr1*k1*V2 + cr2*y2))/
     -   (cr1*k1**2*V1*V2 + cr3*k1*V2*y1 + cr2*k1*V1*y2 + cr4*y1*y2)**2
      gR(2,2)=(k1*V2*y1*((cr1 - cr2)*k1*V1 + (cr3 - cr4)*y1))/
     -   (cr1*k1**2*V1*V2 + cr3*k1*V2*y1 + cr2*k1*V1*y2 + cr4*y1*y2)**2
      gR(3,2)=(k2*V2*(k2*V3 + y3)*(cr5*k2*V3 + cr7*y3))/
     -   (cr5*k2**2*V2*V3 + cr6*k2*V3*y2 + cr7*k2*V2*y3 + cr8*y2*y3)**2
      gR(3,3)=(k2*V3*y2*((cr5 - cr7)*k2*V2 + (cr6 - cr8)*y2))/
     -   (cr5*k2**2*V2*V3 + cr6*k2*V3*y2 + cr7*k2*V2*y3 + cr8*y2*y3)**2
      gR(4,2)=k3/(V2*1.0d+15)
      gR(5,3)=(k2*V3*(k2*V4 + y4)*(cr9*k2*V4 + cr10*y4))/
     -   (cr9*k2**2*V3*V4 + cr11*k2*V4*y3 + cr10*k2*V3*y4 +
     - cr12*y3*y4)**2
      gR(5,4)=(k2*V4*y3*((-cr10 + cr9)*k2*V3 + (cr11 - cr12)*y3))/
     -   (cr9*k2**2*V3*V4 + cr11*k2*V4*y3 + cr10*k2*V3*y4 +
     - cr12*y3*y4)**2
      gR(6,4)=(k4*V4*(k4*V5 + y5)*(cr13*k4*V5 + cr15*y5))/
     -   (cr13*k4**2*V4*V5 + cr14*k4*V5*y4 + cr15*k4*V4*y5 + 
     -   cr16*y4*y5)**2
      gR(6,5)=(k4*V5*y4*((cr13 - cr15)*k4*V4 + (cr14 - cr16)*y4))/
     -   (cr13*k4**2*V4*V5 + cr14*k4*V5*y4 + cr15*k4*V4*y5 + 
     -   cr16*y4*y5)**2    
      gR(7,3)=(k5*V3*Va)/(k5*V3 + y3)**2
      gR(8,3)=fac*(k5*V3*Va)/(12.0d0*(k5*V3 + y3)**2)
      gR(9,5)=(k5*V5*Vn)/(k5*V5 + y5)**2
      gR(10,5)=(k5*V5*Vn)/(12.0d0*(k5*V5 + y5)**2)
      gR(11,6)=(c2*k6*(c2*(cr18 - cr20) + (cr17 - 
     - cr19)*k6)*V1)/(k6*(c2*cr18 + cr17*k6)*V1 + (c2*cr20 +
     - cr19*k6)*y6)**2
      gR(12,6)=(k6*V1*(k6*V2 + y7)*(cr17*k6*V2 + cr18*y7))/
     -   (k6*V2*(cr17*k6*V1 + cr19*y6) + (cr18*k6*V1 + 
     - cr20*y6)*y7)**2
      gR(12,7)=(k6*V2*y6*((cr17 - cr18)*k6*V1 + (cr19 - 
     - cr20)*y6))/(k6*V2*(cr17*k6*V1 + cr19*y6) + 
     - (cr18*k6*V1 + cr20*y6)*y7)**2
      gR(13,7)=(k7*V2*(k7*V3 + y8)*(cr21*k7*V3 + cr23*y8))/
     -   (k7*V3*(cr21*k7*V2 + cr22*y7) + (cr23*k7*V2 + 
     - cr24*y7)*y8)**2
      gR(13,8)=(k7*V3*y7*((cr21 - cr23)*k7*V2 + (cr22 -
     - cr24)*y7))/(k7*V3*(cr21*k7*V2 + cr22*y7) + 
     - (cr23*k7*V2 + cr24*y7)*y8)**2
      gR(14,7)=k3/(V2*1.0d+15) 
      gR(15,8)=(k7*V3*(k7*V4 + y9)*(cr25*k7*V4 + cr26*y9))/
     -   (k7*V4*(cr25*k7*V3 + cr27*y8) + (cr26*k7*V3 +
     - cr28*y8)*y9)**2
      gR(15,9)=(k7*V4*y8*((cr25 - cr26)*k7*V3 + (cr27 -
     - cr28)*y8))/(k7*V4*(cr25*k7*V3 + cr27*y8) +
     - (cr26*k7*V3 + cr28*y8)*y9)**2
      gR(16,9)=(k8*V4*(k8*V5 + y10)*(cr29*k8*V5 + cr31*y10))/
     -   (k8*V4*(cr29*k8*V5 + cr31*y10) + (cr30*k8*V5 +
     - cr32*y10)*y9)**2
      gR(16,10)=(k8*V5*y9*((cr29 - cr31)*k8*V4 + (cr30 -
     - cr32)*y9))/(k8*V4*(cr29*k8*V5 + cr31*y10) +
     - (cr30*k8*V5 + cr32*y10)*y9)**2
      gR(17,8)=(k9*V3*Val)/(k9*V3 + y8)**2
      gR(18,10)=(k9*V5*Vnl)/(k9*V5 + y10)**2
      gR(19,1)=(k1*(c1 + k1)*(c1*cr2 + cr1*k1)*V1)/
     -   (c1*cr2*k1*V1 + cr1*k1**2*V1 + c1*cr4*y1 + cr3*k1*y1)**2
      gR(20,1)=(k1*V1*y2*((cr1 - cr3)*k1*V2 + (cr2 - cr4)*y2))/
     -   (cr1*k1**2*V1*V2 + cr3*k1*V2*y1 + cr2*k1*V1*y2 + cr4*y1*y2)**2
      gR(20,2)=(k1*V2*(k1*V1 + y1)*(cr1*k1*V1 + cr3*y1))/
     -   (cr1*k1**2*V1*V2 + cr3*k1*V2*y1 + cr2*k1*V1*y2 + cr4*y1*y2)**2
      gR(21,2)=(k2*V2*y3*((cr5 - cr6)*k2*V3 + (cr7 - cr8)*y3))/
     -   (cr5*k2**2*V2*V3 + cr6*k2*V3*y2 + cr7*k2*V2*y3 + cr8*y2*y3)**2
      gR(21,3)=(k2*V3*(k2*V2 + y2)*(cr5*k2*V2 + cr6*y2))/
     -   (cr5*k2**2*V2*V3 + cr6*k2*V3*y2 + cr7*k2*V2*y3 + cr8*y2*y3)**2
      gR(22,4)=k3/(V4*1.0d+15)
      gR(23,3)=(k2*V3*y4*((-cr11 + cr9)*k2*V4 + (cr10 - cr12)*y4))/
     -   (cr9*k2**2*V3*V4 + cr11*k2*V4*y3 + cr10*k2*V3*y4 +
     - cr12*y3*y4)**2
      gR(23,4)=(k2*V4*(k2*V3 + y3)*(cr9*k2*V3 + cr11*y3))/
     -   (cr9*k2**2*V3*V4 + cr11*k2*V4*y3 + cr10*k2*V3*y4 +
     - cr12*y3*y4)**2
      gR(24,4)=(k4*V4*y5*((cr13 - cr14)*k4*V5 + (cr15 - cr16)*y5))/
     -   (cr13*k4**2*V4*V5 + cr14*k4*V5*y4 + cr15*k4*V4*y5 + 
     -   cr16*y4*y5)**2
      gR(24,5)=(k4*V5*(k4*V4 + y4)*(cr13*k4*V4 + cr14*y4))/
     -   (cr13*k4**2*V4*V5 + cr14*k4*V5*y4 + cr15*k4*V4*y5 +
     -   cr16*y4*y5)**2
      gR(25,6)=(k6*(c2 + k6)*(c2*cr18 + cr17*k6)*V1)/
     -   (k6*(c2*cr18 + cr17*k6)*V1 + (c2*cr20 +
     - cr19*k6)*y6)**2
      gR(26,6)=(k6*V1*y7*((cr17 - cr19)*k6*V2 + (cr18 -
     - cr20)*y7))/(k6*V2*(cr17*k6*V1 + cr19*y6) +
     - (cr18*k6*V1 + cr20*y6)*y7)**2
      gR(26,7)=(k6*V2*(k6*V1 + y6)*(cr17*k6*V1 + cr19*y6))/
     -   (k6*V2*(cr17*k6*V1 + cr19*y6) + (cr18*k6*V1 +
     - cr20*y6)*y7)**2
      gR(27,7)=(k7*V2*y8*((cr21 - cr22)*k7*V3 + (cr23 -
     - cr24)*y8))/(k7*V3*(cr21*k7*V2 + cr22*y7) +
     - (cr23*k7*V2 + cr24*y7)*y8)**2
      gR(27,8)=(k7*V3*(k7*V2 + y7)*(cr21*k7*V2 + cr22*y7))/
     -   (k7*V3*(cr21*k7*V2 + cr22*y7) + (cr23*k7*V2 +
     - cr24*y7)*y8)**2
      gR(28,9)=k3/(V4*1.0d+15)
      gR(29,8)=(k7*V3*y9*((cr25 - cr27)*k7*V4 + (cr26 -
     - cr28)*y9))/(k7*V4*(cr25*k7*V3 + cr27*y8) +
     - (cr26*k7*V3 + cr28*y8)*y9)**2
      gR(29,9)=(k7*V4*(k7*V3 + y8)*(cr25*k7*V3 + cr27*y8))/
     -   (k7*V4*(cr25*k7*V3 + cr27*y8) + (cr26*k7*V3 +
     - cr28*y8)*y9)**2
      gR(30,9)=(k8*V4*y10*((cr29 - cr30)*k8*V5 + (cr31 -
     - cr32)*y10))/(k8*V4*(cr29*k8*V5 + cr31*y10) +
     - (cr30*k8*V5 + cr32*y10)*y9)**2
      gR(30,10)=(k8*V5*(k8*V4 + y9)*(cr29*k8*V4 + cr30*y9))/
     -   (k8*V4*(cr29*k8*V5 + cr31*y10) + (cr30*k8*V5 +
     - cr32*y10)*y9)**2


      return
      end
cpd-----------------------------------------------
cpd
cpd------------------------------------------------------------------
cpd   SUBROUT used to calulate grad(R) analytically from mathematica
      subroutine GLCRatio(T,N,k,YY,Ras,Rne)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: Ras,Rne
      double precision :: RR(k)
c      double precision :: y1,y2,y3,y4,y5,y6,y7,y8,y9,y10
      integer :: i,j

      
      call rates(n,t,yy,k,RR)
      Ras=(RR(3)-RR(21)-RR(5)+RR(23))/(RR(3)-RR(21)-RR(5)+
     - RR(23)+RR(6)-RR(24))
      Rne=1.0d0-Ras
      
      return
      end