module TimeProp
  use prec
  use couplings
  use hamil

  implicit none

  contains 

  subroutine Propagation(ks, inp)
    implicit none
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp

    integer      :: tion, tele
    integer      :: i, j
    real(kind=q) :: edt
    real(kind=q) :: start, fin

    edt = inp%POTIM / inp%NELM
    ! write(*,*) inp%POTIM, inp%NELM, edt

    ! the OUTER loop
    do tion = 1, inp%NAMDTIME - 1
      ks%pop_a(:,tion) = CONJG(ks%psi_c) * ks%psi_c
      ks%norm(tion)    = SUM(ks%pop_a(:,tion))
      ks%psi_a(:,tion) = ks%psi_c

      
      ! write(*,*) tion, ks%norm(tion), ks%psi_c(:)                     ! check the norm of the state
      ! call cpu_time(start)
      do tele = 1, inp%NELM-1                                           ! the INNER loop
        call make_hamil(tion, tele, ks, inp)                            ! construct hamiltonian matrix
        call hamil_act(ks)                                              ! apply hamiltonian to state vector
        if (tion == 1 .AND. tele == 1) then
          ! write(*,*) ((ks%ham_c(i,j), j=1, ks%ndim), i=1, ks%ndim)
                                                                        ! This is the very first step of the time propagation
                                                                        ! use first order difference
                                                                        ! [c,n,p] means current, next, previous respectively
          ks%psi_n = ks%psi_c - imgUnit * edt * ks%hpsi / hbar
        else
          ks%psi_n = ks%psi_p - 2 * imgUnit * edt * ks%hpsi / hbar      ! use second order difference
        end if
        ks%psi_p = ks%psi_c
        ks%psi_c = ks%psi_n
      end do                                                            ! end of the INNER loop
                                                                        ! call cpu_time(fin)
                                                                        ! write(*,*) "T_ion ", tion, fin - start
    end do                                                              ! end of the OUTER loop
  end subroutine
end module
