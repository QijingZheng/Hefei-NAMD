module hamil
  use prec
  use fileio
  use couplings
  use constants
  implicit none

  type TDKS
    integer :: ndim
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:) :: psi_c
    !!!!Old Integration
    !!!!complex(kind=q), allocatable, dimension(:) :: psi_p
    !!!!complex(kind=q), allocatable, dimension(:) :: psi_n
    complex(kind=q), allocatable, dimension(:,:) :: psi_a
    ! the result of hamiltonian acting on a vector
    complex(kind=q), allocatable, dimension(:) :: hpsi
    ! population
    real(kind=q), allocatable, dimension(:,:) :: pop_a
    real(kind=q), allocatable, dimension(:) :: norm 
    real(kind=q) :: norm_c
    !!!!real(kind=q) :: norm_c, norm_n
    !!!!real(kind=q), allocatable, dimension(:) :: pop_c, pop_n

    complex(kind=q), allocatable, dimension(:,:) :: ham_c
    !!!!complex(kind=q), allocatable, dimension(:,:) :: ham_dt
    ! complex(kind=q), allocatable, dimension(:,:) :: ham_p
    ! complex(kind=q), allocatable, dimension(:,:) :: ham_n

    ! KS eigenvalues
    real(kind=q), allocatable, dimension(:,:) :: eigKs
    ! Non-adiabatic couplings
    complex(kind=q), allocatable, dimension(:,:,:) :: NAcoup

    !! surface hopping related
    !! Bkm = REAL(DCONJG(Akm) * Ckm)
    real(kind=q), allocatable, dimension(:) :: Bkm
    real(kind=q), allocatable, dimension(:,:) :: sh_pops
    real(kind=q), allocatable, dimension(:,:) :: sh_prop

    !! decoherence induced surface hopping
    !!real(kind=q), allocatable, dimension(:,:) :: dish_pops
    !!real(kind=q), allocatable, dimension(:,:) :: recom_pops
    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.

  end type

  contains

  subroutine initTDKS(ks, inp, olap)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(overlap), intent(in)  :: olap
    type(namdInfo), intent(in) :: inp

    integer :: i, j, N

    ! memory allocation
    ks%ndim = inp%NBASIS
    N = inp%NBASIS

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      !!!!allocate(ks%psi_p(N))
      !!!!allocate(ks%psi_n(N))
      allocate(ks%hpsi(N))
      allocate(ks%psi_a(N, 0:inp%NAMDTIME))
      allocate(ks%pop_a(N, 0:inp%NAMDTIME))
      allocate(ks%norm(inp%NAMDTIME))
      !!!!allocate(ks%pop_c(N))
      !!!!allocate(ks%pop_n(N))

      allocate(ks%ham_c(N,N))
      !!!!allocate(ks%ham_dt(N,N))

      !allocate(ks%eigKs(N, inp%NAMDTIME))
      allocate(ks%eigKs(N, inp%NSW))
      !allocate(ks%NAcoup(N, N, inp%NAMDTIME))
      allocate(ks%NAcoup(N, N, 0:inp%NSW-1))

      allocate(ks%sh_pops(N, 0:inp%NAMDTIME))
      allocate(ks%sh_prop(N, 0:inp%NAMDTIME))
      allocate(ks%Bkm(N))
      !!allocate(ks%dish_pops(N, inp%RTIME))
      !!allocate(ks%recom_pops(N,inp%RTIME))
      !!!! allocate(ks%ham_p(N,N))
      !!!! allocate(ks%ham_n(N,N))
      ks%LALLO = .TRUE.
    end if

    ! cero = (0, 0)
    ! uno = (1, 0)
    ks%psi_c = cero
    !!!! ks%psi_p = cero
    !!!! ks%psi_n = cero
    !!!! ks%ham_c = cero
    !!!! ks%ham_p = cero
    !!!! ks%ham_n = cero
    ks%psi_c(inp%INIBAND - inp%BMIN + 1) = uno

    !Using NAMDTIME in NAC loading
 
    do i=1, inp%NAMDTIME
    
    ! We don't need all the information, only a section of it
       ks%eigKs(:,i) = olap%Eig(:, inp%NAMDTINI + i - 1)
    ! Divide by 2 * POTIM here, because we didn't do this in the calculation
    ! of couplings
       ks%NAcoup(:,:,i) = olap%Dij(:,:, inp%NAMDTINI + i - 1) / (2*inp%POTIM)
    end do

    ! new hamil construnction needs nac at t-dt, t, t+dt
    ! so let nac(0-dt) = nac(0)
    if (inp%NAMDTINI == 1) then
       ks%NAcoup(:,:,0) = olap%Dij(:,:,inp%NAMDTINI) / (2*inp%POTIM)
    else   
       ks%NAcoup(:,:,0) = olap%Dij(:,:,inp%NAMDTINI - 1 ) / (2*inp%POTIM)
    end if

    !In DISH, to replicate NAC, we load all NACs.
    !ks%eigKs = olap%Eig
    !ks%NAcoup = olap%Dij / (2*inp%POTIM)
  end subroutine

  ! constructing the hamiltonian by replicating NAC
  subroutine make_hamil_rtime(tion, ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    !real(kind=q) :: dt
    integer :: RTIME,XTIME
    integer :: i
    RTIME = MOD(tion+inp%NAMDTINI-1,inp%NSW-1)
    XTIME = RTIME + 1
    if (RTIME == 0) then
      RTIME = inp%NSW-1
      XTIME = inp%NSW 
    end if
    !XTIME is only used for calculating the average energy now

    ! the hamiltonian contains two parts, which are obtained by interpolation
    ! method between two ionic tims step

  
    ! The non-adiabatic coupling part
    ks%ham_c(:,:) = ks%NAcoup(:,:,RTIME) 
    ! ks%ham_dt(:,:) =(ks%NAcoup(:,:,XTIME) - ks%NAcoup(:,:,RTIME))  / REAL(inp%NELM,q)

    ! multiply by -i * hbar
    ks%ham_c = -imgUnit * hbar * ks%ham_c
    !ks%ham_dt = -imgUnit * hbar * ks%ham_dt

    ! the energy eigenvalue part
    do i=1, ks%ndim
      !Hii(t+0.5dt)
      ks%ham_c(i,i) = 0.5_q*(ks%eigKs(i,RTIME)+ks%eigKs(i,XTIME))
      !ks%ham_c(i,i) = ks%eigKs(i,RTIME) 
      ! ks%ham_dt(i,i) =(ks%eigKs(i,XTIME) - ks%eigKs(i,tion))  / REAL(inp%NELM,q)
    end do

    !if(tion> 967) then
    !    write(967,*) "eig",ks%eigKs(i,tion)
    !end if
  end subroutine

  ! constructing the hamiltonian
  subroutine make_hamil(tion,  ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion

    integer :: i

    ! the hamiltonian contains two parts, which are obtained by interpolation
    ! method between two ionic tims step

    ! The non-adiabatic coupling part
    ks%ham_c(:,:) = ks%NAcoup(:,:,tion)
    !               (ks%NAcoup(:,:,tion+1) - ks%NAcoup(:,:,tion)) * TELE / inp%NELM

    ! multiply by -i * hbar
    ks%ham_c = -imgUnit * hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      !Hii(t+0.5dt)
      ks%ham_c(i,i) = 0.5_q * (ks%eigKs(i,tion) + ks%eigKs(i,tion+1))
      !ks%ham_c(i,i) = ks%eigKs(i,tion)
      !                     (ks%eigKs(i,tion+1) - ks%eigKs(i,tion)) * TELE / inp%NELM
    end do
  end subroutine

  subroutine make_hamil2(tion, TELE,  ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion, TELE

    integer :: i

    ! the hamiltonian contains two parts, which are obtained by interpolation
    ! method between two ionic tims step

    ! The non-adiabatic coupling part
    ks%ham_c(:,:) = ks%NAcoup(:,:,tion) + (ks%NAcoup(:,:,tion+1) - ks%NAcoup(:,:,tion)) * TELE / inp%NELM

    ! multiply by -i * hbar
    ks%ham_c = -imgUnit * hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      !Hii(t+0.5dt)
      !ks%ham_c(i,i) = 0.5_q * (ks%eigKs(i,tion) + ks%eigKs(i,tion+1))
      ks%ham_c(i,i) = ks%eigKs(i,tion) +  (ks%eigKs(i,tion+1) - ks%eigKs(i,tion)) * TELE / inp%NELM
    end do

  end subroutine

  subroutine make_hamil3(tion, TELE,  ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion, TELE

    integer :: i

    ! the hamiltonian contains two parts, which are obtained by interpolation
    ! method between two ionic tims step

    ! The non-adiabatic coupling part
    ! non-adiabatic coupling is actually at (t0+t1)/2
    if (TELE <= (inp%NELM / 2)) then
       ks%ham_c(:,:) = ks%NAcoup(:,:,tion - 1) + (ks%NAcoup(:,:,tion) - ks%NAcoup(:,:,tion-1)) * (TELE+inp%NELM/2 - 0.5_q) / inp%NELM
    else 
       ks%ham_c(:,:) = ks%NAcoup(:,:,tion ) + (ks%NAcoup(:,:,tion+1) - ks%NAcoup(:,:,tion)) * (TELE-inp%NELM/2 - 0.5_q ) / inp%NELM
    end if

    ! multiply by -i * hbar
    ks%ham_c = -imgUnit * hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    ! slighty different, but it won't be a big deal
    do i=1, ks%ndim
      !Hii(t+0.5dt)
      !ks%ham_c(i,i) = 0.5_q * (ks%eigKs(i,tion) + ks%eigKs(i,tion+1))
      ks%ham_c(i,i) = ks%eigKs(i,tion) +  (ks%eigKs(i,tion+1) - ks%eigKs(i,tion)) * (TELE - 0.5_q) / inp%NELM
    end do

  end subroutine



  ! Acting the hamiltonian on the state vector
  subroutine hamil_act(ks)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer :: i, j, N
    !complex(kind=q),dimension(ks%ndim) :: tmp

    N = ks%ndim
    ks%hpsi = cero
    do i=1, N
      do j=1, N
        ks%hpsi(j) =ks%hpsi(j)+ks%ham_c(j,i) * ks%psi_c(i)
      end do
    end do
    !ks%hpsi= matmul(ks%psi_c, ks%ham_c)
  end subroutine

end module
