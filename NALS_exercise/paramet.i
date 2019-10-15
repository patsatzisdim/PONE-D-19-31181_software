	  DOUBLE PRECISION c1,c2,V1,V2,V3,V4,V5
	  DOUBLE PRECISION k1,k2,k3,k4,k5,k6,k7,k8,k9
	  DOUBLE PRECISION cr1,cr2,cr3,cr4,cr5,cr6,cr7,cr8,cr9,cr10,cr11
	  DOUBLE PRECISION cr12,cr13,cr14,cr15,cr16,cr17,cr18,cr19,cr20
	  DOUBLE PRECISION cr21,cr22,cr23,cr24,cr25,cr26,cr27,cr28,cr29
	  DOUBLE PRECISION cr30,cr31,cr32
	  DOUBLE PRECISION GCF,OCF,Va,Vn,Val,Vnl,GRB,ORB
        DOUBLE PRECISION fac,FM


c---------------------------------
c     a parameter to control the flow j8 (R8)
c---------------------------------	  
      parameter (fac=1.0d0)        
c---------------------------------
c     a parameter to change the Rates for evaluating diagnostics
c---------------------------------          
      parameter (FM=1.0d0)
c---------------------------------
c     Syrums and H+
c---------------------------------
      parameter (c1=5.5D0)       !GLC 
      parameter (c2=2.0D+1)       !LAC
c---------------------------------
c     Volumes
c---------------------------------
      parameter (V1=1.2852D-15)   !endothelium
      parameter (V2=1.02D-15)     !basal lamina
      parameter (V3=19.72D-15)    !astrocyte
      parameter (V4=12.24D-15)    !interstitium
      parameter (V5=33.66D-15)    !neuron
c---------------------------------
c     Diffusion Constant
c---------------------------------
c1      parameter (k3=6.0D0)       !kapp in paper
      parameter (k3=1.0D0)       !kapp in supple    !Supple
c---------------------------------
c     Endothelial GLC transp
c---------------------------------
      parameter (cr2=13.0D+13/0.4D0)       !Roi
      parameter (cr3=13.0D+13/0.4D0)       !Rio
      parameter (cr4=13.0D+13/0.4D0)       !Ree
      parameter (cr1=cr2+cr3-cr4)    !Roo
      parameter (k1=8.0D0)          !K
c---------------------------------
c     Astrocyte basal lamina GLC transp    astrocyte to interstitium*12.5
c
c     Change agian for trials!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (decrease of rate=> increase of the trnasport coefficients)
c                                                                *1.5
c---------------------------------
      parameter (cr6=6.25D+15/0.4D0)       !Roi  
      parameter (cr7=5.0D+15/0.4D0)       !Rio
      parameter (cr8=5.0D+15/0.4D0)       !Ree
      parameter (cr5=cr6+cr7-cr8)          !Roo
      parameter (k2=8.0D0)                 !Ka
c---------------------------------
c     Astrocyte interstitium GLC transp
c---------------------------------
      parameter (cr10=5.0D+14/0.4D0)       !Roi
      parameter (cr11=3.8D+14/0.4D0)       !Rio
      parameter (cr12=3.8D+14/0.4D0)       !Ree
      parameter (cr9=cr10+cr11-cr12)       !Roo
c---------------------------------
c     Neuron GLC transp
c---------------------------------
      parameter (cr14=4.4D+13/0.4D0)         !Roi
      parameter (cr15=3.2D+13/0.4D0)         !Rio
      parameter (cr16=3.2D+13/0.4D0)         !Ree
      parameter (cr13=cr14+cr15-cr16)  !Roo
      parameter (k4=4.0D0)             !Kn
c---------------------------------
c     Endothelial lactate transp
c---------------------------------
      parameter (cr18=2.8D+15)         !Roi
      parameter (cr19=2.0D+15)         !Rio
      parameter (cr20=2.0D+15)         !Ree
      parameter (cr17=cr18+cr19-cr20)  !Roo
      parameter (k6=8.0D0)             !eKl  
c---------------------------------
c     Astrocyte basal lamina lactate transp     astrocyte to interstitium*12.5
c
c     Change agian for trials!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (decrease of rate=> increase of the trnasport coefficients)
c                                                                *1.5
c
c     Back to normal
c---------------------------------
      parameter (cr22=8.25D+14)       !Roi
      parameter (cr23=4.125D+14)      !Rio
      parameter (cr24=1.25D+14)       !Ree
      parameter (cr21=cr22+cr23-cr24)     !Roo   
      parameter (k7=5.0d0)                !aKl     
c---------------------------------
c     Astrocyte interstitium lactate transp
c---------------------------------
      parameter (cr26=6.6D+13)          !Roi
      parameter (cr27=3.3D+13)          !Rio
      parameter (cr28=1.0D+13)          !Ree
      parameter (cr25=cr26+cr27-cr28)   !Roo    
c---------------------------------
c     Neuron lactate transp
c---------------------------------
      parameter (cr30=2.0D+14)          !Roi
      parameter (cr31=1.0D+13)          !Rio
      parameter (cr32=1.0D+13)          !Ree
      parameter (cr29=cr30+cr31-cr32)   !Roo 
      parameter (k8=0.7D0)              !nKl    
c---------------------------------
c     Metabolic constants
c---------------------------------
      parameter (k5=4.5D-2)                            !hexokinase MM const
      parameter (k9=2.0D0)                             !lactate MM const
      parameter (GRB=0.0845D0)                         !scaling factor
      parameter (GCF=GRB+0.0d0)
      parameter (Va=GCF/cr10)                 !astrocyte Vmax hexokinase       !Yes
      parameter (Vn=GRB*0.9d0/(2.95D0*cr14))    !neuron Vmax hexokinase          !No change during simulation
      parameter (ORB=0.184D0)                          !scaling factor
      parameter (OCF=ORB+0.0d0)
      parameter (Val=OCF*0.24D0/cr26)         !astrocyte Vmax lactate          !Yes
      parameter (Vnl=1.22D0*OCF/cr30)        !neuron Vmax lactate             !Yes          

