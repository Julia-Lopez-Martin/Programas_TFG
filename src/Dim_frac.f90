Program Frac_dim_v2
    !gfortran -O3 Frac_dim_v2.f90 -o frac.exe && frac.exe
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    ! ------------------PARAMETROS DE AJUSTE -------------------------
    !d_min y d_max en principio se calculan mas adelante, solo hay que decidir Num_dat
    integer, parameter :: Num_dat = 10
    integer, parameter :: modelo_frac = 1
    integer, parameter :: k_max = 1000   
    real(dp):: d_min 
    real(dp):: d_max   
    ! --------------------------------------------

    integer :: P, dim, npar, nfiles
    real(dp) :: dim_fract, c, mu, nu
    real(dp), dimension(:,:), allocatable :: X
    real(dp) :: d(Num_dat)
    integer :: N(Num_dat)
    integer, allocatable :: kpar(:)
    real(dp), allocatable :: mupar(:), nupar(:)
    integer :: i, k, kk, ios
    logical :: existe
    character(len=256) :: fichero, nombre, linea

    dim = 2
    ! Si la carpeta dimension no existe, se crea
    open(unit=9, file="dimension/prueba.tmp", status="replace", iostat=ios)
    if (ios /= 0) then
        call execute_command_line("mkdir dimension")
    else
        close(unit=9, status="delete")
    end if

    ! ---- leer mu y nu de parametros.dat (si existe) ----
    npar = 0
    open(unit=3, file="trayectorias/parametros.dat", status="old", action="read", iostat=ios)
    if (ios == 0) then
        read(3,'(A)',iostat=ios) linea
        do
            read(3,*,iostat=ios) kk, mu, nu, d_min, d_max
            if (ios /= 0) exit
            npar = npar + 1
        end do
        rewind(3)
        read(3,'(A)') linea
        allocate(kpar(npar), mupar(npar), nupar(npar))
        do i = 1, npar
            read(3,*) kpar(i), mupar(i), nupar(i)
        end do
        close(3)
    else
        allocate(kpar(0), mupar(0), nupar(0))
    end if

    open(unit=2, file="dimension/dimensiones_fractales.dat", status="replace", iostat=ios)
    if (ios /= 0) then
        write(*,*) "No se puede crear la carpeta/archivo en dimension/"
        stop
    end if
    write(2,'(A)') "# num   mu   nu   dimension_fractal"

    nfiles = 0
    do k = 1, k_max

        write(fichero,'(A,I0,A)') "trayectorias/trayectoria_", k, ".dat"
        inquire(file=trim(fichero), exist=existe)
        if (.not. existe) cycle

        mu = 0.0_dp; nu = 0.0_dp
        do i = 1, npar
            if (kpar(i) == k) then
                mu = mupar(i); nu = nupar(i)
                exit
            end if
        end do

        open(unit=10, file=trim(fichero), status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write(*,*) "No se puede abrir ", trim(fichero)
            cycle
        end if

        P = 0
        do
            read(10,*,iostat=ios)
            if (ios /= 0) exit
            P = P + 1
        end do
        rewind(10)

        allocate(X(P,dim))
        do i = 1, P
            read(10,*) X(i,1), X(i,2)
        end do
        close(10)
     
        call calc_d(X,d_min,d_max,Num_dat,d,N,P,dim,dim_fract,c)

        ! N frente a d de esta trayectoria
        write(nombre,'(A,I0,A)') "dimension/dimension_", k, ".dat"
        open(unit=1, file=trim(nombre), status="replace")
        write(1,'(A)') "# d   N"
        do i = 1, Num_dat
            write(1,*) d(i), N(i)
        end do
        close(1)

        write(2,*) k, mu, nu, dim_fract
        write(*,*) "Trayectoria", k, " dimension fractal:", dim_fract
        nfiles = nfiles + 1

        deallocate(X)
    end do

    close(2)

    if (nfiles == 0) then
        write(*,*) "No se ha encontrado ninguna trayectoria en trayectorias/ (ejecuta el programa desde la carpeta que la contiene)"
    else
        write(*,*) "Trayectorias analizadas:", nfiles
    end if

    contains

    subroutine Calc_d(X,d_min,d_max,Num_dat,d,N,P,dim,dim_fract,c)
        Implicit none
        integer, intent(in) :: P, dim
        real(dp), intent(in)  :: X(P,dim)
        real(dp), intent(in) :: d_min,d_max
        integer, intent(in) :: Num_dat

        real(dp), intent(out) :: d(Num_dat)
        integer, intent(out) :: N(Num_dat)
        real(dp),intent(out) :: dim_fract,c

        real(dp) :: ld(Num_dat),lN(Num_dat), sumx,sumy,sumxy,sumxx
        integer :: i

        call Calc_N(X,d_min,d_max,Num_dat,d,N,P,dim)

        sumx=0.0_dp
        sumy=0.0_dp
        sumxx=0.0_dp
        sumxy=0.0_dp

        do i=1,Num_dat
            if(N(i)/=0) then
                lN(i)=log(1.0_dp*N(i))
                ld(i)=log(1.0_dp/d(i))

                sumy=sumy+lN(i)
                sumx=sumx+ld(i)
                sumxx=sumxx+ld(i)**2
                sumxy=sumxy+ld(i)*lN(i)
            else
                print*, "Algo  ha salido mal porque ha salido N=0"
            end if
        end do

        dim_fract=(Num_dat*sumxy-sumx*sumy)*1.0_dp/(Num_dat*sumxx-sumx*sumx)
        c=(sumy-dim_fract*sumx)*(1.0_dp/Num_dat)

    end subroutine

    subroutine Calc_N(X,d_min,d_max,Num_dat,d,N,P,dim)
        implicit none
        integer, intent(in) :: P, dim
        real(dp), intent(in)  :: X(P,dim)
        real(dp), intent(in) :: d_min,d_max
        integer, intent(in) :: Num_dat

        real(dp), intent(out) :: d(Num_dat)
        integer, intent(out) :: N(Num_dat)

        integer :: i
        real(dp) :: d_aux

        d_aux= log(d_max/d_min)/(Num_dat-1)

        do i=1,Num_dat
            d(i) = d_min*exp(d_aux*(i-1))
            call calc_Num(X,d(i),N(i),P,dim)
        end do

    end subroutine

    subroutine calc_Num(X, d, N, P, dim)
        implicit none
        integer, intent(in)  :: P, dim
        real(dp), intent(in) :: X(P,dim)
        real(dp), intent(in) :: d
        integer, intent(out) :: N

        integer(8), allocatable :: key(:)
        integer(8) :: gx, gy
        real(dp) :: inv_d
        integer :: i

        allocate(key(P))
        inv_d = 1.0_dp / d

        do i = 1, P
            gx = floor(X(i,1) * inv_d, kind=8)
            gy = floor(X(i,2) * inv_d, kind=8)

            key(i) = (gx + 2_8**30) * 2_8**31 + (gy + 2_8**30)
        end do

        call sort_int8(key, P)

        N = 1
        do i = 2, P
            if (key(i) /= key(i-1)) N = N + 1
        end do

        deallocate(key)
    end subroutine

    recursive subroutine sort_int8(a, n)
        integer(8), intent(inout) :: a(:)
        integer, intent(in) :: n
        integer(8) :: pivot, tmp
        integer :: i, j

        if (n < 2) return
        pivot = a(n/2 + 1)
        i = 1
        j = n
        do
            do while (a(i) < pivot); i = i + 1; end do
            do while (a(j) > pivot); j = j - 1; end do
            if (i >= j) exit
            tmp = a(i); a(i) = a(j); a(j) = tmp
            i = i + 1
            j = j - 1
        end do
        call sort_int8(a(1:i-1), i-1)
        call sort_int8(a(j+1:n), n-j)
    end subroutine

end program