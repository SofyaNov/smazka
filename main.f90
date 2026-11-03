program smazka
    implicit none
    integer :: NI, NJ, max_iter, N_sec, order_acc, section
    integer :: IO_IN, IO_RES, IO_OUT, IO_F
    integer :: i, j
    real(8) :: R_in, eps, L, delta_r, h1, P0, dr, dphi
    real(8) :: F
    real(8), allocatable :: r(:,:), phi(:,:), X(:,:), Y(:,:)
    real(8), allocatable :: h(:,:), p(:,:), p_a(:,:)

    IO_IN=10
    call readinput(IO_IN,"input/Input.txt",NI,NJ,max_iter, N_sec, order_acc,R_in, eps,L,delta_r, h1,P0,dr,dphi,section)
    allocate(r(1:NI,1:NJ), phi(1:NI,1:NJ), X(1:NI,1:NJ), Y(1:NI,1:NJ))
    allocate(h(1:NI,1:NJ), p(1:NI,1:NJ), p_a(1:NI,1:NJ))

    call setka(NI,NJ,dr,dphi,r,R_in,phi,X,Y)
    p = 0.0d0

    IO_RES=11
    call Gauss_Seidel(NI, NJ, max_iter,section, order_acc, N_sec, IO_RES,dr, dphi,r,phi,R_in, delta_r, h1, L, EPS,h,p,F)

    call analitika(NI,NJ,N_sec,L,h1,h,r,p_a)

    call force(Ni,Nj,r,phi,p,F)
    IO_F=13
    open(IO_F, FILE="force.plt", status='unknown', action='write', position='append')
        !write(IO_f,'(I0,1X,I0,1X,I0,1X,ES24.16,1X,ES24.16)') NI,NJ,N_sec,delta_r,F
        write(IO_f,'(I0,1X,I0,1X,I0,1X,G0,1X,G0)') NI,NJ,N_sec,delta_r,F*N_sec
    close(IO_F)

    IO_OUT=12
    call Output_Fields(IO_OUT,"output/fields.plt", NI, NJ, X, Y,r, R_in, phi, H,P,P_a)

    

    deallocate(r, phi, X, Y, h, p,p_a)
end program smazka









subroutine Output_Fields(IO,filename, NI, NJ, X, Y,r, R_in, phi, H,P,P_a)
    implicit none
    integer, intent(in) :: IO, NI, NJ
    real(8), intent(in) :: R_in
    CHARACTER(*), intent(in) :: filename
    real(8), dimension(1:NI,1:NJ), intent(in) :: X, Y,phi, H,P,P_a
    real(8), dimension(1:NI,1:NJ), intent(inout) :: r
    r(:,:)=r(:,:)+R_in
    open(IO, FILE=filename, status='replace', action='write')
        Write(IO,*) 'VARIABLES = "X", "Y", "r", "phi", "H", "P", "P_a"'
        Write(IO,*) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK'
        Write(IO,'(100E25.16)') X(1:NI,1:NJ)
        Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        Write(IO,'(100E25.16)') r(1:NI,1:NJ)
        Write(IO,'(100E25.16)') phi(1:NI,1:NJ)
        Write(IO,'(100E25.16)') H(1:NI,1:NJ)
        Write(IO,'(100E25.16)') P(1:NI,1:NJ)
        Write(IO,'(100E25.16)') P_a(1:NI,1:NJ)
    close(IO)
end subroutine Output_Fields


subroutine readinput(IO,filename,NI,NJ,max_iter, N_sec, order_acc,R_in, eps,L,delta_r, h1,P0,dr,dphi,section)
    implicit none
    INTEGER, intent(in) :: IO
    CHARACTER(*), intent(in) :: filename
    INTEGER, intent(out) :: NI,NJ,max_iter, N_sec, order_acc, section
    REAL(8), intent(out) :: R_in, eps,L,delta_r, h1,P0, dr,dphi
    real(8) :: pi
    pi=4.0d0 * atan(1.0d0)
    open(IO, FILE=filename, status='old', action='read')
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) R_in
        read(IO,*) eps
        read(IO,*) max_iter
        read(IO,*) L
        read(IO,*) delta_r
        read(IO,*) h1
        read(IO,*) N_sec
        read(IO,*) P0
        read(IO,*) order_acc
    close(IO)

    write(*,*) 'R_in  = ', R_in
    write(*,*) 'NI    = ', NI
    write(*,*) 'NJ    = ', NJ
    write(*,*) 'eps   = ', eps
    write(*,*) 'max_iter = ', max_iter
    write(*,*) 'N_sec = ', N_sec
    write(*,*) 'order_acc = ', order_acc

    dr=1.0/(NI-1)
    !dr=(1.0-R_in)/(NI-1)
    !NI=nint((1.0-R_in)/dr)+2
    pi = 4.0d0 * atan(1.0d0)
    dphi = 2.0d0 * (pi/N_sec) / (NJ-1)
    section = NINT(delta_r/dr+1)
    if (abs((section-1)*dr - delta_r) > 1.0d-4) then
        write(*,*) 'delta_r не кратна dr; фактическая ширина =', (section-1)*dr
    end if
END subroutine readinput


real(8) function hf(R_in,delta_r,h1,r,phi,N_sec)
    implicit none
    integer, intent(in) :: N_sec
    real(8), intent(in) :: R_in,delta_r,h1,r,phi
    real(8) :: pi
    pi = 4.0d0 * atan(1.0d0)
    if (r>delta_r+1e-7 .and. r<1.0-delta_r-1e-7) then
        hf=h1+N_sec/2.0d0/pi*(1.0d0-h1)*phi
    else
        hf=1.0d0
    end if
end function hf


subroutine setka(NI,NJ,dr,dphi,r,R_in,phi,X,Y)
    implicit none
    integer,intent(in) :: NI,NJ
    real(8), intent(in) :: dr, dphi, R_in
    real(8), dimension(1:NI,1:NJ), intent(out) :: r,phi,X,Y
    integer :: i,j

    do i = NI,1,-1
        r(i,:) = 1.0 - (NI-i) * dr
    end do
    do j = 1, NJ
        phi(:,j) = (j-1) * dphi
    end do
    X(:,:) = (r(:,:)+R_in) * cos(phi(:,:))
    Y(:,:) = (r(:,:)+R_in) * sin(phi(:,:))
end subroutine setka


subroutine force(Ni,Nj,r,phi,p,F)
implicit none
    integer :: i,j
    integer, intent(in) :: Ni,Nj
    real(8) :: dr_square,dphi
    real(8), dimension(1:Ni,1:Nj) :: r,phi
    real(8), dimension(1:Ni,1:Nj) :: p
    real(8), intent(out) :: F
    F = 0.0d0
    DO I=1,Ni-1
        DO J=1,Nj-1
            dr_square=r(i+1,j+1)*r(i+1,j+1)-r(i,j)*r(i,j)
            dphi=phi(i+1,j+1)-phi(i,j)
            F=F+dr_square*dphi*0.125d0*(P(I,J)+P(I,J+1)+P(I+1,J)+P(I+1,J+1))
        END DO
    END DO
end subroutine force


subroutine boundary_conditions(NI, NJ, h, p, section, order_acc)
    implicit none
    integer, intent(in) :: NI, NJ, section, order_acc
    real(8), dimension(1:NI,1:NJ), intent(in) ::  h
    real(8), dimension(1:NI,1:NJ), intent(inout) :: p
    real(8), dimension(1:NJ) :: znamenatel, chislitel1, chislitel2
    integer :: i,j

    if (order_acc==1) then
        chislitel1=(h(section+1,:))**3*p(section+1,:)+(h(section-1,:))**3*p(section-1,:)
        znamenatel=(h(section-1,:))**3+(h(section+1,:))**3
        p(section,:)=chislitel1/znamenatel
        chislitel1=(h(NI-section+2,:))**3*p(NI-section+2,:)+(h(NI-section,:))**3*p(NI-section,:)
        znamenatel=(h(NI-section,:))**3+(h(NI-section+2,:))**3
        p(NI-section+1,:)=chislitel1/znamenatel
    else if (order_acc==2) then
        chislitel1=4.0d0*(h(section-1,:))**3*p(section-1,:)+4.0d0*(h(section+1,:))**3*p(section+1,:)
        chislitel2=(h(section-1,:))**3*p(section-2,:)+(h(section+1,:))**3*p(section+2,:)
        znamenatel=3.0d0*((h(section-1,:))**3+(h(section+1,:))**3)
        p(section,:)=(chislitel1-chislitel2)/znamenatel
        chislitel1=4.0d0*(h(NI-section,:))**3*p(NI-section,:)+4.0d0*(h(NI-section+2,:))**3*p(NI-section+2,:)
        chislitel2=(h(NI-section,:))**3*p(NI-section-1,:)+(h(NI-section+2,:))**3*p(NI-section+3,:)
        znamenatel=3.0d0*((h(NI-section,:))**3+(h(NI-section+2,:))**3)
        p(NI-section+1,:)=(chislitel1-chislitel2)/znamenatel
    end if

    p(1,:)  = 0.0d0
    p(NI,:) = 0.0d0
    p(:,1)  = 0.0d0
    p(:,NJ) = 0.0d0

end subroutine boundary_conditions


subroutine analitika(NI,NJ,N_sec,L,h1,h,r,p_a)
    implicit none
    integer, intent(in) :: N_sec, NI,NJ
    real(8), intent(in) :: L,h1
    real(8), dimension(1:NI,1:NJ), intent(in) :: h,r
    real(8), dimension(1:NI,1:NJ), intent(out) :: p_a
    real(8) :: pi
    pi = 4.0d0 * atan(1.0d0)
    p_a=2.0d0*pi*L*r**2/N_sec*(h1-h)*(h-1)/(h1**2-1)/h**2
end subroutine analitika










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

    OPEN(IO_RES,FILE='output/residuals.plt')
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
        if (mod(nb,100)==0) then 
            WRITE(*,*) nb, max_err/max_err1, N_sec*F
        end if
        WRITE(IO_RES,*) nb, max_err/max_err1, N_sec*F
        IF ((max_err/max_err1 <= EPS) .OR. (nb >= max_iter)) THEN
            EXIT
        ENDIF
        p_old=p
    end do
    close(IO_RES)
end subroutine Gauss_Seidel
