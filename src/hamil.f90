module hamil
  use prec
  use fileio
  use couplings
  use constants
  implicit none

  type TDKS 
    integer                                      :: ndim        ! # of states selected
                                                                ! ndim == BMAX - BMIN == NBASIS
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:)   :: psi_c
    complex(kind=q), allocatable, dimension(:)   :: psi_p
    complex(kind=q), allocatable, dimension(:)   :: psi_n
    complex(kind=q), allocatable, dimension(:,:) :: psi_a       ! All the states selected, psi_all
                                                                ! psi_a(NBASIS, NAMDTIME)
    complex(kind=q), allocatable, dimension(:)   :: hpsi        ! the result of hamiltonian acting on a vector
    ! population
    real(kind=q), allocatable, dimension(:,:)    :: pop_a       ! Population of each state
                                                                ! pop_a(NBASIS, NAMDTIME)
    real(kind=q), allocatable, dimension(:)      :: norm        ! normalization factor
                                                                ! Psi = norm * psi_a ?

    complex(kind=q), allocatable, dimension(:,:) :: ham_c       ! Hamiltonian matrix
    ! complex(kind=q), allocatable, dimension(:,:) :: ham_p
    ! complex(kind=q), allocatable, dimension(:,:) :: ham_n
    
    ! KS eigenvalues
    real(kind=q), allocatable, dimension(:,:)    :: eigKs       ! EigVal of KS orbit extracted from WAVECAR
    real(kind=q), allocatable, dimension(:,:,:)  :: NAcoup      ! Non-adiabatic couplings

    ! surface hopping related

    ! Bkm = REAL(CONJG(Akm) * Ckm)
    real(kind=q), allocatable, dimension(:)      :: Bkm
    real(kind=q), allocatable, dimension(:,:)    :: sh_pops
    real(kind=q), allocatable, dimension(:,:)    :: sh_prop

    ! whether the memory has been allocated
    ! Avoid fuzzy meaning, changed var name: LALLO -> LALLOCATED
    logical                                      :: LALLOCATED = .FALSE.

  end type

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initializes time-dependent kohn-sham wavefunction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initTDKS(ks, inp, olap)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(overlap), intent(in)  :: olap
    type(namdInfo), intent(in) :: inp

    integer :: i, j, N

    ! memory allocation
    ks%ndim = inp%NBASIS      ! NBASIS == BMAX - BMIN == # of bands considered
    N       = inp%NBASIS      ! N == NBASIS

    if (.NOT. ks%LALLOCATED) then  ! Don't allocated memory ??
      allocate(ks%psi_c(N))
      allocate(ks%psi_p(N))
      allocate(ks%psi_n(N))
      allocate(ks%hpsi(N))
      allocate(ks%psi_a(N, inp%NAMDTIME))
      allocate(ks%pop_a(N, inp%NAMDTIME))
      allocate(ks%norm(inp%NAMDTIME))

      allocate(ks%ham_c(N,N))

      allocate(ks%eigKs(N, inp%NAMDTIME))
      allocate(ks%NAcoup(N,N, inp%NAMDTIME))

      allocate(ks%sh_pops(N, inp%NAMDTIME))
      allocate(ks%sh_prop(N, inp%NAMDTIME))
      allocate(ks%Bkm(N))
      ! allocate(ks%ham_p(N,N))
      ! allocate(ks%ham_n(N,N))
      ks%LALLOCATED = .TRUE.
    end if

    ! cero = (0, 0)
    ! uno = (1, 0)
    ks%psi_c = cero
    ks%psi_p = cero
    ks%psi_n = cero
    ! ks%ham_c = cero
    ! ks%ham_p = cero
    ! ks%ham_n = cero
    ks%psi_c(inp%INIBAND - inp%BMIN + 1) = uno

    do i=1, inp%NAMDTIME
      ! We don't need all the information, only a section of it
      ks%eigKs(:,i) = olap%Eig(:, inp%NAMDTINI + i - 1)
      ! Divide by 2 * POTIM here, because we didn't do this in the calculation
      ! of couplings
      ks%NAcoup(:,:,i) = olap%Dij(:,:, inp%NAMDTINI + i - 1) / (2*inp%POTIM)
    end do

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! constructing the hamiltonian
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_hamil(TION, TELE, ks, inp)
    implicit none

    type(TDKS), intent(inout)  :: ks                ! Kohn-Sham orbit
    type(namdInfo), intent(in) :: inp               ! Input options
    integer, intent(in)        :: TION, TELE        ! TION <- ionic time step
                                                    ! TELE <- electronic time index

    integer :: i

    ! the hamiltonian contains two parts, which are obtained by interpolation
    ! method between two ionic tims step

    ! The non-adiabatic coupling part
    ks%ham_c(:,:) = ks%NAcoup(:,:,TION) + &
                   (ks%NAcoup(:,:,TION+1) - ks%NAcoup(:,:,TION)) * TELE / inp%NELM

    ! multiply by -i * hbar
    ks%ham_c = -imgUnit * hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = ks%eigKs(i,TION) + &                                        ! Diagonal part is the eigvals from WAVECAR
                     (ks%eigKs(i,TION+1) - ks%eigKs(i,TION)) * TELE / inp%NELM
    end do
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Applying the hamiltonian on the state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine hamil_act(ks)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer                   :: i, j, N
    complex(kind=q)           :: tmp

    N = ks%ndim
    do i=1, N
      tmp = cero                                    ! cero <- complex zero === (0, 0)
      do j=1, N
        tmp = tmp + ks%ham_c(i,j) * ks%psi_c(j)     ! hpsi = H * |psi_c>
      end do
      ks%hpsi(i) = tmp                              ! why not use ks%hpsi = matmul(ks&ham_c, ks%psi_c) ?
    end do

  end subroutine

end module
