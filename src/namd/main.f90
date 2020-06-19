Program main
  use prec
  use lattice
  use wavecar
  use couplings
  use hamil
  use shop
  use fileio
  use TimeProp

  implicit none

  type(namdInfo) :: inp
  type(TDKS) :: ks
  type(overlap) :: olap, olap_sec

  real(kind=q) :: start, fin
  integer :: ns

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, get user inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getUserInp(inp)
  ! call printUserInp(inp)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Secondly, get couplings
  ! In the very first run, the following subroutine will calculate the
  ! NA-couplings from WAVECARs and then write it to a binary file called
  ! COUPCAR.  From the second run on, the subroutine will just read the
  ! NA-couplings from the file. However, for a general NAMD run, the file is way
  ! too huge, the solution is to write only the information we need to another
  ! plain text file. If such files exist (set LCPTXT = .TRUE. in the inp), then
  ! we may skip the huge binary file and read the plain text file instead. This
  ! is done in the subroutine 'initTDKS'.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call TDCoupIJ(trim(inp%rundir), olap, olap_sec, inp)
  ! write(*,*) "T_coup: ", fin - start
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ns=1, inp%NSAMPLE
    inp%NAMDTINI = inp%NAMDTINI_A(ns)
    inp%INIBAND  = inp%INIBAND_A(ns)
    call printUserInp(inp)
    ! initiate KS matrix
    call cpu_time(start)
    call initTDKS(ks, inp, olap_sec)
    ! Time propagation
    call Propagation(ks, inp)
    ! Run surface hopping
    if (inp%LSHP) then
      call runSH(ks, inp)
      call printSH(ks, inp)
    end if
    call cpu_time(fin)
    write(*,'(A, F8.2)') "CPU Time [s]:", fin - start
  end do

end Program
