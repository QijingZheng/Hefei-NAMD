module dish
  use prec
  use fileio
  use hamil
  use TimeProp
  implicit none

  contains

  ! initialize the random seed from the system clock
  ! code from: http://fortranwiki.org/fortran/show/random_seed
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

  subroutine calcdect(ks, inp, DECOTIME, COEFFISQ)
    implicit none
    type(TDKS), intent(in) :: ks
    type(namdInfo), intent(in) :: inp

    real(kind=q), dimension(inp%NBASIS), intent(inout) :: DECOTIME
    real(kind=q), dimension(inp%NBASIS), intent(inout) :: COEFFISQ
    integer :: i, j
    real(kind=q) :: tmp

    do j=1, inp%NBASIS
      COEFFISQ(j) = DCONJG(ks%psi_c(j)) * ks%psi_c(j)
    end do

    do i=1, inp%NBASIS
      tmp = 0
      do j=1, inp%NBASIS
        tmp = tmp + inp%DEPHMATR(j, i) * COEFFISQ(j)   ! DEPHMATR(j,i)=DEPHMATR(i,j)
      end do
      DECOTIME(i) = 1.0_q / tmp
    end do
    !write(90,*),DECOTIME
  end subroutine

  subroutine whichtodec(tion, inp, DECOTIME, which, decmoment,shuffle)
    implicit none
    type(namdInfo), intent(in) :: inp
    real(kind=q), intent(in) :: DECOTIME(inp%NBASIS)
    integer, intent(in) :: tion
    integer, intent(inout) :: which
    real(kind=q), dimension(inp%NBASIS), intent(inout) :: decmoment
    integer, intent(inout),dimension(inp%NBASIS) :: shuffle
    integer :: tmp
    integer :: i, j,n
    integer :: tm(inp%NBASIS), decstate(inp%NBASIS)
    !! the colomn of decstate is the sum number of decoherence states (n)
    !! the element of decstate is the serial number of decoherence states (which)
    real(kind=q) :: r

    n = 0
    decstate = 0
    which = 0


    do j=inp%NBASIS,1,-1
      call random_number(r)
      i=INT(inp%NBASIS*r)+1
      tmp=shuffle(i)
      shuffle(i)=shuffle(j)
      shuffle(j)=tmp
    end do


    do j=1,inp%NBASIS
      i=shuffle(j)
      !tm(i) = tion - decmoment(i)
      if ( DECOTIME(i) <= decmoment(i) ) then
        which = i
        !decmoment(which) = tion    !! update the decoherence moment
        decmoment(which) = 0.0_q    !! update the decoherence moment

        exit
      end if
    end do
  end subroutine

  subroutine projector(ks, inp, COEFFISQ, which, tion, cstat,iend,fgend)
    implicit none
    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    real(kind=q), intent(in) :: COEFFISQ(inp%NBASIS)
    integer, intent(in) :: which, tion,iend
    integer, intent(inout) :: cstat,fgend

    integer :: i, j
    real(kind=q) :: r, dE, kbT, popmax,lower,upper
    real(kind=q) :: popBoltz(inp%NBASIS)
    integer :: RTIME 

    RTIME = MOD(TION+inp%NAMDTINI-1,inp%NSW-1)
    if (RTIME == 0) then
      RTIME = inp%NSW-1
    end if

    kbT = inp%TEMP * BOLKEV
    call random_number(r)

    !! hopping probability with Boltzmann factor
    popBoltz(which) = DCONJG(ks%psi_c(which)) * ks%psi_c(which)
    !if (Mod(tion,100)==0) then
    !  write(89,*) which,tion,popBoltz(which)
    !end if
    if (inp%LHOLE) then
      dE = ((ks%eigKs(cstat,RTIME) + ks%eigKs(cstat,RTIME+1)) - &
              (ks%eigKs(which,RTIME) + ks%eigKs(which, RTIME+1))) /2.0_q
      if ( dE > 0.0_q ) then
        popBoltz(which) = popBoltz(which) * exp(-dE / kbT)
      end if
    else
      dE = ((ks%eigKs(which,RTIME) + ks%eigKs(which,RTIME+1)) - &
              (ks%eigKs(cstat,RTIME) + ks%eigKs(cstat, RTIME+1))) /2.0_q
      if ( dE > 0.0_q ) then
        popBoltz(which) = popBoltz(which) * exp(-dE / kbT)
      end if
    end if
    !write(*,*) which,popBoltz(which)
    !! project in/out
    if (r <= popBoltz(which)) then
      do i=1, inp%NBASIS
        ks%psi_c(i) = cero
      end do
      ks%psi_c(which) = uno

      if (fgend==0 .AND. which==iend) then
        ks%recom_pops(cstat,tion + 1:inp%RTIME)=ks%recom_pops(cstat, tion + 1:inp%RTIME)+1.0_q
        fgend=-1
      end if

      cstat = which

    else
    !! project out
      ks%psi_c(which) = cero
      ks%norm_c=SUM(DCONJG(ks%psi_c)*ks%psi_c)
      !write(83,*) tion,which,ks%psi_c,ks%norm_c
      do i=1, inp%NBASIS
        ks%psi_c(i) = ks%psi_c(i) / DSQRT(ks%norm_c)
      end do
      !ks%psi_c(which) = cero
      !if (cstat==which) then
      !    call random_number(r)
      !    do i=1,inp%NBASIS
      !        if (i==1) then
      !            lower = 0.0_q
      !            upper = DCONJG(ks%psi_c(i)) * ks%psi_c(i)
      !        else
      !            lower = upper
      !            upper = upper+DCONJG(ks%psi_c(i)) * ks%psi_c(i)

      !        end if

      !        if(lower <= r .AND. r< upper) then
      !            if (fgend==0 .AND. i==iend) then

      !                ks%recom_pops(cstat,tion + 1:inp%RTIME)=ks%recom_pops(cstat, tion + 1:inp%RTIME)+1.0_q
      !                fgend=-1
      !            end if
      !            cstat=i
      !            exit
      !        end if

      !    end do

      !end if
    end if

  end subroutine

  subroutine runDISH(ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer :: i, tion,k,j
    integer :: istat, cstat, which,iend,fgend

    real(kind=q) :: COEFFISQ(inp%NBASIS)
    real(kind=q) :: DECOTIME(inp%NBASIS)
    integer, dimension(inp%NBASIS) :: decmoment
    integer, dimension(inp%NBASIS) :: shuffle

    ks%dish_pops = 0.0_q
    ks%recom_pops = 0.0_q
    istat = inp%INIBAND - inp%BMIN + 1
    if (inp%LHOLE) then
      iend=inp%BMAX-inp%BMIN+1
    else
      iend=1
    end if


    ! initialize the random seed for ramdom number production
    call init_random_seed()
    !$OMP PARALLEL DO SHARED(ks,inp,istat), PRIVATE(I,cstat),REDUCTION(+:ks%dish_pops)
    !MPI is under development
    do i=1, inp%NTRAJ
      fgend=0
      call mpirundish(ks,inp,istat,cstat,iend,fgend)
    end do
    !$OMP END PARALLEL DO
    ks%dish_pops = ks%dish_pops / inp%NTRAJ
    ks%dish_pops(istat, 1) = 1.0_q

  end subroutine

  subroutine mpirundish(ks,inp,istat,cstat,iend,fgend)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(inout) :: fgend
    integer, intent(in) :: iend
    integer :: i, tion,k,j
    integer, intent(inout) :: istat, cstat
    integer :: which
    real(kind=q), dimension(inp%NBASIS) :: decmoment
    integer, dimension(inp%NBASIS) :: shuffle

    integer :: start,fin 
    real(kind=q) :: COEFFISQ(inp%NBASIS)
    real(kind=q) :: DECOTIME(inp%NBASIS)
    ! At the first step, current state always equal initial state
    ks%psi_c = cero
    ks%psi_c(istat) = uno
    cstat = istat
    decmoment = 0.0_q

    do j=1,inp%NBASIS
      shuffle(j)=j
    end do

    !do tion=1, inp%NAMDTIME - 1
    do tion=1, inp%RTIME - 1
      !call system_clock(start)
      call PropagationT(ks, inp, tion)
      !call system_clock(fin)
      !write(*,'(A, I)') "CPU Time in Propagation [s]:", fin - start
      !start=fin      
      call calcdect(ks, inp, DECOTIME, COEFFISQ)
      !call system_clock(fin)
      !write(*,'(A, I)') "CPU Time in calcdect [s]:", fin - start
      !start=fin      

      call whichtodec(tion, inp, DECOTIME, which, decmoment,shuffle)
      !call system_clock(fin)
      !write(*,'(A, I)') "CPU Time in whichtodec [s]:", fin - start
      !start=fin      
        
      decmoment=decmoment+inp%POTIM
      if ( which > 0 ) then
        call projector(ks, inp, COEFFISQ, which, tion, cstat,iend,fgend)
      end if
      ks%dish_pops(cstat, tion + 1) = ks%dish_pops(cstat, tion + 1) + 1.0_q

    end do

  end subroutine

  subroutine printDISH(ks, inp)
    implicit none
    type(TDKS), intent(in) :: ks
    type(namdInfo), intent(in) :: inp

    integer :: i, j, tion, ierr, io
    character(len=48) :: buf
    character(len=48) :: out_fmt
    integer :: RTTIME
    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f9.4))" )' )  ks%ndim

    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=28, file='RECOMB.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    !! open(unit=25, file='PSICT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "SHPROP file I/O error!"
      stop
    end if

    !! do io = 24 !!, 25
    io = 24
    write(io,'(A,A12,A3,I5)') '#', 'BMIN',     ' = ', inp%BMIN
    write(io,'(A,A12,A3,I5)') '#', 'BMAX',     ' = ', inp%BMAX
    write(io,'(A,A12,A3,I5)') '#', 'INIBAND',  ' = ', inp%INIBAND
    write(io,'(A,A12,A3,I5)') '#', 'NBANDS',   ' = ', inp%NBANDS

    write(io,'(A,A12,A3,I5)')   '#', 'NSW',    ' = ', inp%NSW
    write(io,'(A,A12,A3,F5.1)') '#', 'POTIM',  ' = ', inp%POTIM
    write(io,'(A,A12,A3,F5.1)') '#', 'TEMP',   ' = ', inp%TEMP

    write(io,'(A,A12,A3,I5)') '#', 'NAMDTINI', ' = ', inp%NAMDTINI
    write(io,'(A,A12,A3,I5)') '#', 'NAMDTIME', ' = ', inp%NAMDTIME
    write(io,'(A,A12,A3,I9)') '#', 'RTIME',    ' = ', inp%RTIME
    write(io,'(A,A12,A3,I5)') '#', 'NTRAJ',    ' = ', inp%NTRAJ
    write(io,'(A,A12,A3,I5)') '#', 'NELM',     ' = ', inp%NELM

    write(io,'(A,A12,A3,L5)') '#', 'LHOLE',    ' = ', inp%LHOLE
    !! write(io,'(A,A12,A3,L5)') '#', 'LSHP',     ' = ', inp%LSHP
    write(io,'(A,A12,A3,L5)') '#', 'LCPTXT',   ' = ', inp%LCPTXT
    write(io,'(A,A12,A3,A)')  '#', 'RUNDIR',   ' = ', TRIM(ADJUSTL(inp%rundir))
    !! end do

    !do tion=1, inp%NAMDTIME
    do tion=1, inp%RTIME
      RTTIME=mod(tion+inp%NAMDTINI-1,inp%NSW-1)
      if (RTTIME==0) then
        RTTIME = inp%NSW-1
      end if


      write(unit=24,fmt=out_fmt) tion * inp%POTIM, SUM(ks%eigKs(:,RTTIME) * ks%dish_pops(:,tion)), &
                            (ks%dish_pops(i,tion), i=1, ks%ndim)
      !! write(unit=25, fmt=*) tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%pop_a(:,tion)), &
                            !! (ks%psi_a(i,tion), i=1, ks%ndim)
                            ! (ks%pop_a(i,tion), i=1, ks%ndim)

      write(unit=28,fmt=out_fmt) tion * inp%POTIM, SUM(ks%eigKs(:,RTTIME) * ks%dish_pops(:,tion)), &
                            (ks%recom_pops(i,tion), i=1, ks%ndim)
    end do

    close(24)
    close(28)
    !! close(25)

  end subroutine

end module
