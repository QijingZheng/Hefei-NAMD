! MODULE for lattice
! 
module lattice
  use prec
  type latt
    real(q) :: SCALE
    ! A: real space basis
    ! B: momentum sapce basis
    real(q) :: A(3,3),B(3,3)
    ! NORM of the basises
    real(q) :: ANORM(3),BNORM(3)
    ! volume of the cell
    real(q) :: OMEGA
  end type
  contains

    ! subroutine for calculating the reciprocal lattice from the direct lattice in
    ! addition the norm of the lattice-vectors and the volume of the basis-cell is
    ! calculated

    subroutine lattic(Mylatt)
      use prec

      implicit none

      type(latt) :: Mylatt
      real(q) :: Omega
      integer :: i,j
      intrinsic SUM

      ! calculate reciprocal space basis from real space one
      CALL EXPRO(Mylatt%B(:,1),Mylatt%A(:,2),Mylatt%A(:,3))
      CALL EXPRO(Mylatt%B(:,2),Mylatt%A(:,3),Mylatt%A(:,1))
      CALL EXPRO(Mylatt%B(:,3),Mylatt%A(:,1),Mylatt%A(:,2))

      ! the volume: B1 \dot A1
      Omega = SUM(Mylatt%A(:,1) * Mylatt%B(:,1))
      ! Omega = Mylatt%B(1,1) * Mylatt%A(1,1) + Mylatt%B(2,1) * Mylatt%A(2,1) &
      ! &       + Mylatt%B(3,1) * Mylatt%A(3,1)

      Mylatt%B = Mylatt%B / Omega
      ! do i=1,3
      !   do j=1,3
      !     Mylatt%B(i,j) = Mylatt%B(i,j) / Omega
      !   end do
      ! end do

      ! magnitude of basis vector
      do i=1,3
        Mylatt%ANORM(i) = SQRT(SUM(Mylatt%A(:,i)*Mylatt%A(:,i)))
        Mylatt%BNORM(i) = SQRT(SUM(Mylatt%B(:,i)*Mylatt%B(:,i)))
      end do
      Mylatt%OMEGA = Omega

      return
    end subroutine

    ! caclulates the x-product of two vectors
    ! H = U1 x U2

    subroutine EXPRO(H,U1,U2)
      use prec
      implicit none
      real(kind=q), intent(in) :: U1, U2
      real(kind=q), intent(inout) :: H
      dimension H(3), U1(3), U2(3)

      H(1)=U1(2)*U2(3)-U1(3)*U2(2)
      H(2)=U1(3)*U2(1)-U1(1)*U2(3)
      H(3)=U1(1)*U2(2)-U1(2)*U2(1)

      return
    end subroutine

end module lattice
