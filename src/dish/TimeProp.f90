module TimeProp
  use prec
  use couplings
  use hamil

  implicit none

  contains 
  
  !Use Trotter formula to integrate the Time-dependtent Schrodinger equation
  !The new scheme provides a robust and efficient integration. The electronic
  !time step (POTIM/NELM) may be increased up to the value of the nuclear time step (Not in all)
  !Always check convergence before using a small NELM 
  !This scheme is proposed by Akimov, A. V., & Prezhdo, O. V. J. Chem. Theory Comput. 2014, 10, 2, 789–804
  
  
  !This scheme has been revised by Dr.Li yunhai (liyunhai1016@hotmail.com)
  !To use this scheme, the offdiagonal elements of Hamiltonian should be real numbers (without
  !the the imaginary unit) 
  
   subroutine PropagationT(ks, inp, tion)
    implicit none
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion

    integer :: tele
    integer :: i, j
    real(kind=q) :: edt
    real(kind=q) :: norm
    

    integer :: jj, kk
    complex(kind=q) :: phi, cos_phi, sin_phi, cjj, ckk
    complex(kind=q), parameter :: miuno = (0.0_q, -1.0_q)
    
    
    



    edt = inp%POTIM / inp%NELM
    ! write(*,*) inp%POTIM, inp%NELM, edt

    ! the OUTER loop
    !!!do tion = 1, inp%NAMDTIME - 1
    !!!  ks%pop_a(:,tion) = CONJG(ks%psi_c) * ks%psi_c
    !!!  ks%norm(tion) = SUM(ks%pop_a(:,tion))
    !!!  ks%psi_a(:,tion) = ks%psi_c

      ! check the norm of the state
      ! write(*,*) tion, ks%norm(tion), ks%psi_c(:)
      ! call cpu_time(start)
      ! the INNER loop
      ! do tele = 1, inp%NELM-1
      do tele = 1, inp%NELM
        ! construct hamiltonian matrix
        !Exact
        call make_hamil_rtime(tion, tele, ks, inp)
        !Convential 
        !call make_hamil2(tion, tele, ks, inp)


        ! propagate the psi_c according to Liouville-Trotter algorithm
        !   exp[i*(L_{ij}+L_i)*dt/2]
        ! = exp(i*L_{ij}*dt/2) * exp(i*L_i*dt/2) * exp(i*L_i*dt/2) * exp(i*L_{ij}*dt/2)
        ! = exp(i*L_{ij}*dt/2) * exp(i*L_i*dt) * exp(i*L_{ij}*dt/2)

        ! First L_{ij} part
        
        !Changed the matrix structure for NAC.
        !Now Ham_c(i,j) is stored as  ks%ham_c(j,i)
        !Weibin

        do jj = 1, inp%NBASIS
          do kk = jj+1, inp%NBASIS
            phi = 0.5_q * edt * miuno * ks%ham_c(kk, jj) / hbar
            cos_phi = cos(phi)
            sin_phi = sin(phi)
            cjj = ks%psi_c(jj)
            ckk = ks%psi_c(kk)
            ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
            ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
          end do
        end do
        
        ! Li part

        !Changed the matrix structure for NAC.
        !Now Ham_c(i,j) is stored as  ks%ham_c(j,i)
        !Weibin


        do jj = 1, inp%NBASIS
          phi = edt * miuno * ks%ham_c(jj, jj) / hbar
          ks%psi_c(jj) = ks%psi_c(jj) * exp(phi)
        end do

        ! Second L_{ij} part
        do jj = inp%NBASIS, 1, -1
          do kk = inp%NBASIS, jj+1, -1
            phi = 0.5_q * edt * miuno * ks%ham_c(kk, jj) / hbar
            cos_phi = cos(phi)
            sin_phi = sin(phi)
            cjj = ks%psi_c(jj)
            ckk = ks%psi_c(kk)
            ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
            ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk

          end do
        end do
      end do
     norm = REAL(SUM(CONJG(ks%psi_c) * ks%psi_c)) 
     if ( norm <= 0.99_q .OR. norm >= 1.01_q)  then
        write(*,*) "Error in Electronic Propagation"
        stop
     end if
 
     !ks%psi_a(:,tion) = ks%psi_c
     !ks%pop_a(:,tion) = ks%pop_a(:,tion)/ks%norm(tion) 


     !write(89,*) tion, conjg(ks%psi_c)*ks%psi_c
     !write(89,*) 'sum', sum(conjg(ks%psi_c)*ks%psi_c)
     ! ! end of the INNER loop
     ! ! call cpu_time(fin)
     ! ! write(*,*) "T_ion ", tion, fin - start
    ! end do
    ! end of the OUTER loop
  end subroutine
 
  


!!  subroutine Propagation(ks, inp, tion)
!!    implicit none
!!    type(TDKS), intent(inout)  :: ks
!!    type(namdInfo), intent(in) :: inp
!!    integer, intent(in) :: tion
!!
!!    integer :: tele
!!    integer :: i, j
!!    real(kind=q) :: edt
!!    real(kind=q) :: start, fin
!!    complex(kind=q) ::ALPHA1,ALPHA2,BETA 
!!    ALPHA1=imgUnit*edt/hbar*(-1.0_q)
!!    ALPHA2=imgUnit*edt/hbar*(-2.0_q)
!!    BETA=(1.0_q,0.0_q)
!!        
!!    edt = inp%POTIM / inp%NELM
!!    ! write(*,*) inp%POTIM, inp%NELM, edt
!!
!!    ! the OUTER loop
!!    !! do tion = 1, inp%NAMDTIME - 1
!!    !! ks%pop_a(:,tion) = DCONJG(ks%psi_c) * ks%psi_c
!!    !! ks%norm(tion) = SUM(ks%pop_a(:,tion))
!!    !! ks%psi_a(:,tion) = ks%psi_c
!!
!!    ! check the norm of the state
!!    ! write(*,*) tion, ks%norm(tion), ks%psi_c(:)
!!    ! call cpu_time(start)
!!    ! the INNER loop
!!    ! do tele = 1, inp%NELM-1
!!    call make_hamil(tion, ks, inp)
!!    !write(98,'(*(E))') tion,real(conjg(ks%psi_c)*ks%psi_c),real(sum(conjg(ks%psi_c)*ks%psi_c))
!!    do tele = 1, inp%NELM-1
!!      ! construct hamiltonian matrix
!!      ks%ham_c=ks%ham_c+ks%ham_dt
!!      !call make_hamil(tion, tele, ks, inp,edt)
!!      ! apply hamiltonian to state vector
!!      !write(*,*) ks%ham_c    
!!      !if (tion > 968) then
!!      !    write(966,*) "psi_c",tion,ks%ham_c,"\n"
!!      !end if
!!        
!!      call hamil_act(ks)
!!      !call ZGEMV(TRANS, M, N, ALPHA, A, LDA, X,INCX, BETA, Y, INCY)
!!      !call ZSYMV('N', inp%NBASIS, inp%NBASIS, ALPHA, A, LDA, X,INCX, BETA, Y, INCY)
!!      if (tion == 1 .AND. tele == 1) then
!!      ! write(*,*) ((ks%ham_c(i,j), j=1, ks%ndim), i=1, ks%ndim)
!!      ! This is the very first step of the time propagation
!!      ! use first order difference
!!      ! [c,n,p] meas current, next, previous respectively
!!      !ks%psi_n=ks%psi_c
!!      !call ZGEMV('N', inp%NBASIS, inp%NBASIS, ALPHA1, ks%ham_c,inp%NBASIS, ks%psi_c,1, BETA, ks%psi_n, 1)
!!        ks%psi_n = ks%psi_c - imgUnit * edt * ks%hpsi / hbar
!!        else
!!      ! use second order difference
!!      !call ZGEMV('N', inp%NBASIS, inp%NBASIS, ALPHA2, ks%ham_c,inp%NBASIS, ks%psi_c,1, BETA, ks%psi_n, 1)
!!        ks%psi_n = ks%psi_p - 2.0_q * imgUnit * edt * ks%hpsi / hbar
!!      end if
!!      !ks%pop_c(:) = DCONJG(ks%psi_c) * ks%psi_c
!!      !ks%norm_c = SUM(ks%pop_c(:))
!!      !ks%psi_p = ks%psi_c / SQRT(ks%norm_c)
!!      !ks%pop_n(:) = CONJG(ks%psi_n) * ks%psi_n
!!      !ks%norm_n = SUM(ks%pop_n(:))
!!      !ks%psi_c = ks%psi_n / SQRT(ks%norm_n)
!!      !if (tion > 968) then
!!      !    write(968,*) "psi_c",tion,tele,ks%psi_c,"\n"
!!      !    write(969,*) "hbar",ks%hpsi,"\n"
!!      !end if
!!      ks%psi_p = ks%psi_c
!!      ks%psi_c = ks%psi_n
!!      !write(92,*) "psi_c",tion,tele,ks%psi_c,"\n"
!!    end do
!!      !write(95,*) "psi_c",tion,tele,conjg(ks%psi_c)*ks%psi_c
!!      ! end of the INNER loop
!!      ! call cpu_time(fin)
!!      ! write(*,*) "T_ion ", tion, fin - start
!!
!!    !! end do
!!    ! end of the OUTER loop
!!  end subroutine

end module
