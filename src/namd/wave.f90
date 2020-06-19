module wavecar
  use prec
  use lattice

  implicit none

  type psi
    integer :: iband
    integer :: ikpts
    integer :: ispin
  end type

  type waveinfo
    ! name of wavecar
    character(len=255) :: WAVECAR
    real(kind=q) :: ENCUT
    type(latt) :: MyLatt
    integer :: ISPIN
    integer :: NKPTS
    integer :: NBANDS
    ! integer :: NGPTAR(3)

    ! precision of wavecar
    integer :: qw = qs
    ! WAVECAR IO unit, default to 13
    integer :: IU

    ! Not needed in coefficients reading
    ! integer, allocatable, dimension(:) :: LPCTX, LPCTY, LPCTZ

    ! number of plane waves for each k-point
    integer, allocatable, dimension(:) :: NPLWS
    ! Maximum plane wave numbers
    integer :: MAXPLWS
    ! k-point vector
    real(kind=q), allocatable, dimension(:,:) :: VKPTS
    real(kind=q), allocatable, dimension(:,:,:) :: BANDS

    ! Not needed in coefficients reading
    ! G-vector index of stored wavefunction
    ! GX[NPLWS, NKPTS]
    ! integer, allocatable, dimension(:,:) :: GX, GY, GZ
  end type

  contains

    subroutine freemem(MySys)
      implicit none
      type(waveinfo), intent(inout) :: MySys

      deallocate(MySys%VKPTS)
      deallocate(MySys%BANDS)
      deallocate(MySys%NPLWS)
      ! deallocate(MySys%GX,MySys%GY,MySys%GZ)
      ! deallocate(MySys%LPCTX,MySys%LPCTY,MySys%LPCTZ)
    end subroutine
    
    subroutine setKet(ket, ib, ik, is)
      implicit none
      type(psi), intent(inout) :: ket
      integer, intent(in) :: ib, ik, is

      ket%iband = ib
      ket%ikpts = ik
      ket%ispin = is
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine to open wavecar
    ! just open "WAVECAR", do not read
    ! associate it to IO unit: IU (12, if IU not specified)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine openwav(wavecar, IU)
      implicit none

      character(len=*), intent(in)   :: wavecar
      integer, intent(in)   :: IU
      integer, parameter :: irecl = 48
      real(kind=q)       :: rdum, rispin, rtag
      integer            :: ierr, idum

      open(unit=IU, file=wavecar, access='direct', form='unformatted', &
           status = 'unknown', recl=irecl, iostat=ierr)
      if(ierr /= 0) then
          write(*,*) "File I/O error with " // wavecar
          stop
      end if

      rdum = 0._q
      read(unit=IU, rec=1, iostat=ierr) rdum, rispin, rtag  ! So, this is the reason why
                                                                ! irecl is 24 at first

      if (NINT(rtag) /= 45200) then
          write(*,*) "WAVECAR not single precision"
      end if

      ! write(*,*) rdum, rispin
      close(unit=IU)

      ! rdum contains the record length
      idum = NINT(rdum)
      if(ierr /= 0 .OR. idum <= 0) then
          write(*,*) "Error reading WAVECAR!"
          stop
      end if
      ! reopen the "WAVECAR" after getting the right RECORD LENGTH
      open(unit=IU, file=wavecar, access='direct', form='unformatted', &
           status = 'unknown', recl=idum, iostat=ierr)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! close unit 12 (usually connected to the  WAVECAR file)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine closewav(IU)
      implicit none
      integer, intent(in) :: IU
      close(IU)
    end subroutine closewav

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! rewind unit 12 (usually connected to the  WAVECAR file)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rewindwav(IU)
      implicit none
      integer, intent(in) :: IU
      rewind(IU)
    end subroutine rewindwav

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read system basic infos from WAVECAR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sysinfo(MySys)
      implicit none

      type(waveinfo), intent(inout) :: MySys

      real(kind=q) :: rdum, rispin, rtag, rnkpts, rnbands, rnplw
      integer :: i, j, k, n, irec
      complex(kind=q), allocatable, dimension(:) :: eig
      real(kind=q), allocatable, dimension(:) :: occ

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! open WAVECAR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call openwav(MySys%WAVECAR, MySys%IU)

      read(unit=MySys%IU, rec=1) rdum, rispin, rtag
      read(unit=MySys%IU, rec=2) rnkpts, rnbands, MySys%ENCUT, ((MySys%MyLatt%A(i,j), i=1,3),j=1,3)

      MySys%ISPIN  = NINT(rispin)
      MySys%NKPTS  = NINT(rnkpts)
      MySys%NBANDS = NINT(rnbands)
      
      if (NINT(rtag) /= 45200) then
          MySys%qw = q
        else
          MySys%qw = qs
      end if

      allocate(MySys%VKPTS(3,MySys%NKPTS))
      allocate(MySys%NPLWS(MySys%NKPTS))
      allocate(MySys%BANDS(MySys%NBANDS,MySys%NKPTS,MySys%ISPIN))
      allocate(eig(MySys%NBANDS))
      allocate(occ(MySys%NBANDS))

      irec = 2
      do i=1, MySys%ISPIN
        do j=1, MySys%NKPTS
          irec = irec + 1
          read(MySys%IU, rec=irec) rnplw, (MySys%VKPTS(k,j), k=1,3), &
                             (eig(n), occ(n), n=1, MySys%NBANDS)  
          
          MySys%NPLWS(j) = NINT(rnplw)
          do n=1, MySys%NBANDS
            MySys%BANDS(n,j,i) = real(eig(n))
          end do
          irec = irec + MySys%NBANDS
        end do
      end do
      MySys%MAXPLWS = MAXVAL(MySys%NPLWS)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! close WAVECAR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! call closewav(MySys%IU)

    end subroutine sysinfo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Load the coefficients of one wavefunction from WAVECAR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine LOADWAVE(cwork, ket, MySys)
      implicit none

      type(waveinfo), intent(in) :: MySys
      type(psi), intent(in) :: ket
      complex(kind=qs), intent(inout) :: cwork(MySys%NPLWS(ket%ikpts))

      integer :: i, irec
      logical :: lopen
      real(kind=q) :: norm

      inquire(file=MySys%WAVECAR, opened=lopen)
      if (.NOT. lopen) then
        write(*,*) "ERROR! WAVECAR not opened..."
        stop
      end if

      irec = 2 + (ket%ispin - 1) * (MySys%NKPTS * (MySys%NBANDS + 1)) + &
                 (ket%ikpts - 1) * (MySys%NBANDS + 1) + (ket%iband + 1)
      
      ! write(*,*) irec, ket%ispin, ket%ikpts, ket%iband
      read(unit=MySys%IU, rec=irec) (cwork(i), i=1,MySys%NPLWS(ket%ikpts))
      ! NORM of the wave function 
      norm = SQRT(SUM(CONJG(cwork) * cwork))
      ! Normalizd the wave function
      cwork = cwork / norm

    end subroutine LOADWAVE

end module wavecar
