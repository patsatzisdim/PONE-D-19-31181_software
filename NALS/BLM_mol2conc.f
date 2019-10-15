cpd   Convert concentrations to molar concentrations
      subroutine mol2conc(n,t,yval,yvalm)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: t, yval(n) ! concentrations in mmol
      double precision, intent(out) :: yvalm(n+2) ! molar concentrations mmol/L
      integer :: i
      double precision :: Vt
      include 'paramet.i'
cpd
cpd   Glucose molar concentrations          
      yvalm(1)=yval(1)/V1    ! endothelium
      yvalm(2)=yval(2)/V2    ! basal lamina
      yvalm(3)=yval(3)/V3    ! astrocyte
      yvalm(4)=yval(4)/V4    ! interstitium
      yvalm(5)=yval(5)/V5    ! neuron
cpd      
cpd   Lactate molar concentrations          
      yvalm(6)=yval(6)/V1    ! endothelium
      yvalm(7)=yval(7)/V2    ! basal lamina
      yvalm(8)=yval(8)/V3    ! astrocyte
      yvalm(9)=yval(9)/V4    ! interstitium
      yvalm(10)=yval(10)/V5    ! neuron
cpd
cpd   total molar concentrations
      yvalm(11)=0.0D0
      yvalm(12)=0.0D0
      do i=1,5
       yvalm(11)=yvalm(11)+yval(i)
       yvalm(12)=yvalm(12)+yval(i+5)
      enddo
      Vt=V1+V2+V3+V4+V5
      yvalm(11)=yvalm(11)/Vt
      yvalm(12)=yvalm(12)/Vt
cpd
      return
      end  