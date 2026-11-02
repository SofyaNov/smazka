program smazka
    implicit none
    integer :: NI, NJ, max_iter, N_sec, order_acc, section
    integer :: IO_IN, IO_RES, IO_OUT
    integer :: i, j
    real(8) :: R_in, eps, L, delta_r, h1, P0, dr, dphi
    real(8) :: F
    real(8), allocatable :: r(:,:), phi(:,:), X(:,:), Y(:,:)
    real(8), allocatable :: h(:,:), p(:,:)

    IO_IN=10
    call readinput(IO_IN,"Input.txt",NI,NJ,max_iter, N_sec, order_acc,R_in, eps,L,delta_r, h1,P0,dr,dphi,section)
    allocate(r(1:NI,1:NJ), phi(1:NI,1:NJ), X(1:NI,1:NJ), Y(1:NI,1:NJ))
    allocate(h(1:NI,1:NJ), p(1:NI,1:NJ))

    call setka(NI,NJ,dr,dphi,r,R_in,phi,X,Y)
    p = 0.0d0

    IO_RES=11
    call Gauss_Seidel(NI, NJ, max_iter,section, order_acc, N_sec, IO_RES,dr, dphi,r,phi,R_in, delta_r, h1, L, EPS,h,p,F)

    IO_OUT=12
    call Output_Fields(IO_OUT, 'fields.plt', NI, NJ, X, Y,r, R_in,phi, P, H)

    deallocate(r, phi, X, Y, h, p)
end program smazka
