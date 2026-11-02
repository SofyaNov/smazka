subroutine Output_Fields(IO,filename, NI, NJ, X, Y,r, R_in, phi, P, H)
    implicit none
    integer, intent(in) :: IO, NI, NJ
    real(8), intent(in) :: R_in
    CHARACTER(*), intent(in) :: filename
    real(8), dimension(1:NI,1:NJ), intent(in) :: X, Y,phi, P, H
    real(8), dimension(1:NI,1:NJ), intent(inout) :: r
    r(:,:)=r(:,:)+R_in
    open(IO, FILE=filename, status='replace', action='write')
        Write(IO,*) 'VARIABLES = "X", "Y", "r", "phi", "P", "H"'
        Write(IO,*) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK'
        Write(IO,'(100E25.16)') X(1:NI,1:NJ)
        Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
        Write(IO,'(100E25.16)') r(1:NI,1:NJ)
        Write(IO,'(100E25.16)') phi(1:NI,1:NJ)
        Write(IO,'(100E25.16)') P(1:NI,1:NJ)
        Write(IO,'(100E25.16)') H(1:NI,1:NJ)
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
