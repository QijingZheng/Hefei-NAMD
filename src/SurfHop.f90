module shop
  use prec
  use fileio
  use hamil
  implicit none

  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize the random seed from the system clock
  ! code from: http://fortranwiki.org/fortran/show/random_seed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_random_seed()
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate `which` state the electron/hole should hop to
  !   `which` is the state index
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine whichToHop(tion, ks, which)
    implicit none

    integer, intent(in)    :: tion              ! ionic time step
    integer, intent(inout) :: which             !
    type(TDKS), intent(in) :: ks

    integer                :: i
    real(kind=q)           :: lower, upper, r

    which = 0
    call random_number(r)

    do i=1, ks%ndim
      if (i == 1) then
        lower = 0
        upper = ks%sh_prop(i,tion)
      else
        lower = upper
        upper = upper + ks%sh_prop(i,tion)
      end if
      if (lower <= r .AND. r < upper) then
        which = i
        exit
      end if
    end do

  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate the hopping probabilities
  ! Gjk = Pjk * Bjk                             Total SH probability imposed on SH action
  ! Pjk = SH probability from FSSH
  ! Bjk = | exp(-deltaE / kBT) if deltaE >  0
  !       |         1          if deltaE <= 0   This part imposes detailed balance to SH action
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calcprop(tion, cstat, ks, inp)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in)        :: tion                          ! ionic time step
    integer, intent(in)        :: cstat                         ! current state index

    integer                    :: i, j
    real(kind=q)               :: Akk
    real(kind=q)               :: dE, kbT

    Akk = CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(cstat, tion)                  ! |ks%psi_a(cstat, tion)|^2
    ! Bkm = REAL(CONJG(Akm) * Ckm)
    ks%Bkm = -2. * REAL(CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(:, tion) * &    ! 
                    ks%NAcoup(cstat, :, tion))                                  ! CPA here

    ks%sh_prop(:,tion) = ks%Bkm / Akk * inp%POTIM

    kbT = inp%TEMP * BOLKEV                                                     ! Boltzmann factor

    if (inp%LHOLE) then
      do i = 1, cstat
        dE = ks%eigKs(cstat, tion) - ks%eigKs(i,tion)
        ks%sh_prop(i,tion) = ks%sh_prop(i,tion) * exp(-dE / kbT)
      end do
    else
      do i=cstat, ks%ndim
        dE = ks%eigKs(i,tion) - ks%eigKs(cstat, tion)
        ks%sh_prop(i,tion) = ks%sh_prop(i,tion) * exp(-dE / kbT)
      end do
    end if

    forall (i=1:ks%ndim, ks%sh_prop(i,tion) < 0) ks%sh_prop(i,tion) = 0
    ! write(*,*) (ks%Bkm(i), i=1, ks%ndim) 
    ! write(*,*) (ks%sh_prop(i, tion), i=1, ks%ndim) 

    !!!!!!
    ! POSSIBLE OPTIMIZATION HERE.
    ! sh_prop is independent of time, could be calculated and stored in advance.
    !!!!!!

  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate surface hopping probabilities
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine runSH(ks, inp)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    integer                    :: i, j, tion                        ! tion  <-  dummy variable, denote time step index
    integer                    :: istat, cstat, which               ! istat  <-  # of states, aka # of bands selected, also be the initial step
                                                                    ! cstat  <-  current stat index
                                                                    ! which  <-  target state that electron/hole will hop to

    ks%sh_pops = 0
    ks%sh_prop = 0
    istat      = inp%INIBAND - inp%BMIN + 1

    ! initialize the random seed for ramdom number production
    call init_random_seed()

    do i=1, inp%NTRAJ
      ! in the first step, current step always equal initial step
      cstat = istat
      do tion=1, inp%NAMDTIME
        call calcprop(tion, cstat, ks, inp)                         ! random stuff not needed, can be calculated in advanced
        call whichToHop(tion, ks, which)                            ! random stuff needed, cannot be predicted
        if (which > 0) cstat = which
        ks%sh_pops(cstat, tion) = ks%sh_pops(cstat, tion) + 1
      end do
    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ
    ! ks%sh_prop = ks%sh_prop / inp%NTRAJ

    ! do tion=1, inp%NAMDTIME
    !   write(*,*) (ks%sh_pops(i,tion), i=1, ks%ndim)
    ! end do
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! need a subroutine here to write the results we need
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine printSH(ks, inp)
    implicit none
    type(TDKS), intent(in) :: ks
    type(namdInfo), intent(in) :: inp

    integer :: i, j, tion, ierr, io
    character(len=48) :: buf

    write(buf, *) inp%NAMDTINI
    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=25, file='PSICT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "SHPROP file I/O error!"
      stop
    end if

    do io = 24, 25
      write(io,'(A,A12,A3,I5)') '#', 'BMIN',     ' = ', inp%BMIN
      write(io,'(A,A12,A3,I5)') '#', 'BMAX',     ' = ', inp%BMAX
      write(io,'(A,A12,A3,I5)') '#', 'INIBAND',  ' = ', inp%INIBAND
      write(io,'(A,A12,A3,I5)') '#', 'NBANDS',   ' = ', inp%NBANDS

      write(io,'(A,A12,A3,I5)')   '#', 'NSW',    ' = ', inp%NSW
      write(io,'(A,A12,A3,F5.1)') '#', 'POTIM',  ' = ', inp%POTIM
      write(io,'(A,A12,A3,F5.1)') '#', 'TEMP',   ' = ', inp%TEMP

      write(io,'(A,A12,A3,I5)') '#', 'NAMDTINI', ' = ', inp%NAMDTINI
      write(io,'(A,A12,A3,I5)') '#', 'NAMDTIME', ' = ', inp%NAMDTIME
      write(io,'(A,A12,A3,I5)') '#', 'NTRAJ',    ' = ', inp%NTRAJ
      write(io,'(A,A12,A3,I5)') '#', 'NELM',     ' = ', inp%NELM

      write(io,'(A,A12,A3,L5)') '#', 'LHOLE',    ' = ', inp%LHOLE
      write(io,'(A,A12,A3,L5)') '#', 'LSHP',     ' = ', inp%LSHP
      write(io,'(A,A12,A3,L5)') '#', 'LCPTXT',   ' = ', inp%LCPTXT
      write(io,'(A,A12,A3,A)')  '#', 'RUNDIR',   ' = ', TRIM(ADJUSTL(inp%rundir))
    end do

    do tion=1, inp%NAMDTIME
      write(unit=24, fmt=*) tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%sh_pops(:,tion)), &
                            (ks%sh_pops(i,tion), i=1, ks%ndim)
      write(unit=25, fmt=*) tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%pop_a(:,tion)), &
                            (ks%psi_a(i,tion), i=1, ks%ndim)
                            ! (ks%pop_a(i,tion), i=1, ks%ndim)
    end do

    close(24)
    close(25)

  end subroutine

end module
