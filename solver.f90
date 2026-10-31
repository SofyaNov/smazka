subroutine Gauss_Seidel(NI, NJ, max_iter,section, order_acc, N_sec, IO_RES,dr, dphi,r,phi,R_in, delta_r, h1, L, EPS,h,p,F)
    implicit none
    integer, intent(in) :: NI, NJ, max_iter,section, order_acc, IO_RES, N_sec
    real(8), intent(in) :: dr, dphi
    real(8), dimension(1:NI,1:NJ), intent(in) :: r,phi
    real(8), intent(in) :: R_in, delta_r, h1, L, EPS
    integer :: i,j, nb
    real(8), dimension(1:NI,1:NJ), intent(inout) :: p
    real(8), dimension(1:NI,1:NJ), intent(out) :: h
    real(8), dimension(1:NI-1,1:NJ) :: r_half
    real(8), dimension(1:NI,1:NJ-1) :: h_half
    real(8) :: A,B,C,D,Q, max_err, max_err1
    real(8), intent(out) :: F
    real(8) :: hf
    real(8), dimension(1:NI,1:NJ) :: p_old

    OPEN(IO_RES,FILE='residuals.plt')
    WRITE(IO_RES,*) 'VARIABLES = "Iterations", "Residuals", "FORCE"'

    do i=1,NI-1
        r_half(i,:)=(r(i,:)+r(i+1,:))/2.0d0
    end do
    do i=1,NI
        do j=1,NJ-1
            h(i,j)=hf(R_in,delta_r,h1,r(i,j),phi(i,j),N_sec)
            h_half(i,j)=hf(R_in,delta_r,h1,r(i,j),(phi(i,j+1)+phi(i,j))/2.0d0,N_sec)
        end do
        h(i,NJ)=hf(R_in,delta_r,h1,r(i,NJ),phi(i,NJ),N_sec)
    end do

    p_old=p

    do nb = 1, max_iter
        max_err = 0.0d0
        do i=2,NI-1
            if (i/=section .and. i/=NI-section+1) then
                do j=2,NJ-1
                    A=(h(i,j))**3/r(i,j)/dr**2*r_half(i,j)
                    B=(h(i,j))**3/r(i,j)/dr**2*r_half(i-1,j)
                    C=(h_half(i,j))**3/r(i,j)**2/dphi**2
                    D=(h_half(i,j-1))**3/r(i,j)**2/dphi**2
                    Q=L*(h(i,j+1)-h(i,j-1))/2.0d0/dphi
                    p(i,j)=(A*p(i+1,j)+B*p(i-1,j)+C*p(i,j+1)+D*p(i,j-1)-Q)/(A+B+C+D)
                end do
            end if
        end do
        call boundary_conditions(NI, NJ, h, p, section, order_acc)

        max_err = MAXVAL(ABS(P_OLD(:,:) - P(:,:)))
        if (nb == 1) THEN
            max_err1=max_err
        end if

        call force(Ni,Nj,r,phi,p,F)
        WRITE(*,*) nb, max_err/max_err1, N_sec*F
        WRITE(IO_RES,*) nb, max_err/max_err1, N_sec*F
        IF ((max_err/max_err1 <= EPS) .OR. (nb >= max_iter)) THEN
            EXIT
        ENDIF
        p_old=p
    end do
    close(IO_RES)
end subroutine Gauss_Seidel

