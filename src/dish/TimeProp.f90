module TimeProp
  use prec
  use couplings
  use hamil

  implicit none

  contains 
  
  !Use Trotter formula to integrate the Time-dependtent Schrodinger equation
  !The new scheme provides a robust and efficient integration. The electronic
  !time step (POTIM/NELM) may be increased up to the value of the nuclear time step
  !Always check convergence before using a small NELM 
  !This scheme is proposed by Akimov, A. V., & Prezhdo, O. V. J. Chem. Theory Comput. 2014, 10, 2, 789–804
  
  subroutine rot1(ks,phi,i,j)
    implicit none
    type(TDKS), intent(inout)  :: ks
    real(kind=q),intent(in):: phi
    integer,intent(in) :: i,j
    complex(kind=q) :: psi_i,psi_j

    psi_i=dcos(phi)*ks%psi_c(i)+dsin(phi)*ks%psi_c(j)
    psi_j=-1.0_q*dsin(phi)*ks%psi_c(i)+dcos(phi)*ks%psi_c(j)
    ks%psi_c(i)=psi_i
    ks%psi_c(j)=psi_j
  end subroutine
  
  subroutine rot2(ks,phi,i,j)
    implicit none
    type(TDKS), intent(inout)  :: ks
    real(kind=q),intent(in):: phi
    integer,intent(in) :: i,j
    complex(kind=q) :: psi_i,psi_j

    psi_i=dcos(phi)*ks%psi_c(i)+imgUnit*dsin(phi)*ks%psi_c(j)
    psi_j=imgUnit*dsin(phi)*ks%psi_c(i)+dcos(phi)*ks%psi_c(j)
    ks%psi_c(i)=psi_i
    ks%psi_c(j)=psi_j
  end subroutine




  subroutine Phase(ks,dt,i)
    implicit none
    type(TDKS), intent(inout)  :: ks
    real(kind=q),intent(in) :: dt 
      
    real(kind=q)  :: phi
    integer,intent(in) :: i

    phi=-dt*real(ks%ham_c(i,i),q)/hbar
    ks%psi_c(i)=(dcos(phi)+imgUnit*dsin(phi))*ks%psi_c(i)
  end subroutine

  subroutine rot(ks,dt,i,j)
    type(TDKS), intent(inout)  :: ks
    real(kind=q),intent(in) :: dt
    integer,intent(in) :: i,j


    real(kind=q) :: phi1,phi2
    
    phi1=0.5_q*dt*DIMAG(ks%ham_c(i,j))/hbar
    phi2=-dt*Real(ks%ham_c(i,j),q)/hbar
  
  
  
    call rot1(ks,phi1,i,j)
    call rot2(ks,phi2,i,j)
    call rot1(ks,phi1,i,j)
  
  
  end subroutine

  subroutine PropagationT(ks,inp,tion)
    implicit none
    type(TDKS), intent(inout)  :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    
    integer :: i,j,tele
    real(kind=q) :: edt

    call make_hamil(tion, ks, inp)
    !Debug only
    !write(98,'(*(E))') tion,real(dconjg(ks%psi_c)*ks%psi_c),real(sum(dconjg(ks%psi_c)*ks%psi_c))
    
    edt = inp%POTIM / inp%NELM
    do tele = 1, inp%NELM
    do i=1,ks%NDIM
      do j= i+1,ks%NDIM
        call rot(ks,0.5_q*edt,i,j)
      end do
    end do

    do i=1,ks%NDIM
      call phase(ks,edt,i)
    end do

    do i=ks%NDIM,1,-1
      do j=ks%NDIM,i+1,-1
        call rot(ks,0.5_q*edt,i,j)
      end do
    end do
    
    !Debug only
    !write(98,*)  "psi",tion,tele,dconjg(ks%psi_c)*ks%psi_c,sum(dconjg(ks%psi_c)*ks%psi_c)
    
    end do   
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
