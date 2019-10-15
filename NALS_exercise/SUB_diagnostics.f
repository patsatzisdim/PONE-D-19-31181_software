cpd------------------------------------------------------------------
cpd   Diagnostics subroutine. Printings are passed out at the main code.
cpd------------------------------------------------------------------
      subroutine Diag(n,k,NoEM,t,yy,alpha,betta,reim,po,tpi,bsr,impi,
     - fip)
      implicit none
      integer, intent(in) :: n,k,NoEM
      integer :: i,j,ii,count,countSP,count1
      double precision, intent(in) :: t, yy(n), reim(n,2)
      double precision, intent(inout) ::  alpha(n,n), betta(n,n)
      double precision, intent(out) :: po(n,n),tpi(n,k),bsr(n,k),
     - fip(n),impi(n,k)
      double precision :: dpn(n), st(n,k), gR(k,n), gradSRk(n,n),
     - tempS(n,1), tempGR(1,n), tempB(1,n), BgSRk(1,n), sum1, 
     - BgSRkA(n,k), tsum1, tsum2, sum2(n), RR(k), FI(4,n),
     - tempSR(n,k), dumBSR(n,k), alphaS(n,n-NoEM), bettaS(n-NoEM,n),
     - asbs(n,n), asbsh(n,k)
      
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
cpd------------------------------
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
cpd------------------------------
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
cpd------------------------------
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
cpd------------------------------      
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
cpd   CSP II through eigencvectors
      call stoic(n,k,st)
      call rates(n,t,yy,k,RR)
      do i=1,n-NoEM
       alphaS(:,i)=alpha(:,i+NoEM)
       bettaS(i,:)=betta(i+NoEM,:)
      enddo
      asbs(:,:)=0.0d0  
      asbsh(:,:)=0.0d0   
cpd------------------------------
      call smult(n,n-NoEM,n,alphaS,bettaS,asbs,n,n-NoEM,n)
      call smult(n,n,k,asbs,st,asbsh,n,n,k)
      do i=1,n
       do j=1,k
        impi(i,j)=asbsh(i,j)*RR(j)
       enddo
       tsum1=0.0d0
       do j=1,k
        tsum1=tsum1+dabs(impi(i,j))+1.0d-20 ! handle with care
       enddo
       if(i.eq.10) then
        if(t.lt.1110.1d0.and.t.gt.111.0d+1) then
         write(*,*) T, impi(i,16), impi(i,30), impi(i,7), impi(i,18),
     -   impi(i,17), impi(i,11), impi(i,9), impi(i,25)
        endif
       endif
       do j=1,k
        impi(i,j)=impi(i,j)/tsum1
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
       	alpha(j,i)=alpha(j,i)*dsign(1.d0,FI(1,i))     ! change the sign of eigenvectors
        betta(j,i)=betta(j,i)*dsign(1.d0,FI(1,i))
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
cpd--------------------------------------------------
c
c      Print the original diagnostics in different files
c
cpd--------------------------------------------------
cpd--------------------------------------------------
      subroutine Print_Diag(n,k,t,yy,po,tpi,api,impi,fip,unNo1,unNo2,
     - unNo3,unNo4,unNo5)
      implicit none
      integer,intent(in) :: n, k, unNo1, unNo2, unNo3, unNo4, unNo5
      double precision, intent(in) :: t, yy(n), po(n,n), api(n,k),
     - tpi(n,k), impi(n,k), fip(n)
      integer :: i

      do i=1,n
       write(unNo1,36) t, 'Spec  ',i, po(i,1:n) 
       write(unNo2,33) t, 'Mode  ',i, tpi(i,1:k)
       write(unNo3,33) t, 'Mode  ',i, api(i,1:k)
       write(unNo4,33) t, 'Spec  ',i, impi(i,1:k)
      enddo
      write(unNo5,30) t, fip(1:n)

   30 format(E18.11,10('    ',E25.16))
   33 format(E18.11,a8,i3,30('    ',E25.16))
   36 format(E18.11,a8,i3,10('    ',E25.16))         

      return
      end  
cpd--------------------------------------------------
c
c      Print the sorted diagnostics in excel format
c
cpd--------------------------------------------------
cpd--------------------------------------------------
      subroutine smart_print(n,k,t,yy,PoS,tpiS,apiS,impiS,PoI,tpiI,
     - apiI,impiI,unNoStart,unNoEnd)
      implicit none
      integer,intent(in) :: n, k, unNoStart, unNoEnd, PoI(n,n), 
     - tpiI(n,k), apiI(n,k), impiI(n,k) 
      double precision, intent(in) :: t, yy(n), PoS(n,n), apiS(n,k),
     - tpiS(n,k), impiS(n,k)
      integer :: i,ii,j

      do ii=unNoStart,unNoEnd         !number of modes
       j=ii-unNoStart+1 
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP pointer'
       do i=1,n
        write(ii,40) t,i,PoI(j,i),PoS(j,i)
       enddo
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP API'
       do i=1,k 
        write(ii,40) t,i,apiI(j,i),apiS(j,i)
       enddo
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP II'       
       do i=1,k
        write(ii,40) t,i,impiI(j,i),impiS(j,i)
       enddo
       write(ii,*) '-----------------------'
       write(ii,*) 'CSP TPI'       
       do i=1,k
        write(ii,40) t,i,tpiI(j,i),tpiS(j,i)
       enddo
      enddo


   40 format(E18.11,i5,i5,F12.3)     

      return
      end 
