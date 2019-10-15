      subroutine sortPRow(n,k,A,B,Ind)
      implicit none
      integer, intent(in) :: n,k   !dimensions of the matrix
      double precision, intent(in) :: A(n,k)    ! matrix to be sorted per row
      double precision, intent(out) :: B(n,k)   ! matrix after sorting 
      integer, intent(out) :: Ind(n,k)          ! matrix for the initial index position
      integer :: i,j,indexIn,ii
      double precision :: max, Row(k)
       
      do i=1,n
       indexIn=0
       do j=1,k
       	Ind(i,j)=indexIn+1
       enddo
      enddo

      do i=1,n                    ! loop over the rows
       do j=1,k
       	Row(j)=A(i,j)
       enddo
       max=Row(1)
       indexIn=1
       do j=1,k
       	do ii=j+1,k
         if(dabs(max).lt.dabs(Row(ii))) then 
          max=Row(ii)
          indexIn=ii
         endif
        enddo
        B(i,j)=max
        Ind(i,j)=indexIn
       enddo
      enddo 


      return
      end

