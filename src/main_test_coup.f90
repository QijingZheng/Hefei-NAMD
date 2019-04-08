Program main
  use prec
  use lattice
  use wavecar
  use couplings

  implicit none

  type(overlap) :: olap
  integer :: i, j, k
  integer :: bandmin, bandmax

  olap%NBANDS = 388
  olap%TSTEPS = 2000
  olap%dt = 1

  ! bandmin = 127
  ! bandmax = 176

  allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS-1))
  allocate(olap%Eig(olap%NBANDS, olap%TSTEPS-1))
  call TDCoupIJ('RUN', olap)

  ! do k=1, 2
  !   do i=bandmin, bandmax
  !     ! Cij = -0.658218 * Cij / 2.
  !     write(*,*) (-0.658218 /2. * olap%Dij(i,j,k), j=bandmin, bandmax)
  !   end do
  !   write(*,*)
  ! end do

end Program
