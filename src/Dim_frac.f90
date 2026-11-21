Program Frac_dim_v2
    !gfortran -O3 Frac_dim_v2.f90 -o frac.exe && frac.exe
    use rutas
    use sort
    
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    ! ------------------PARAMETROS DE AJUSTE -------------------------
    !d_min y d_max en principio se calculan mas adelante, solo hay que decidir Num_dat y modelo
    integer, parameter :: Num_dat = 10

    !modelo 2, numero de puntos desde los que empezar a contar
    integer :: Pr
    integer ::

    !Numero máximo de datos, no cambiar
    integer, parameter :: k_max = 1000   
    real(dp):: d_min 
    real(dp):: d_max   
    ! --------------------------------------------

    integer :: P, dim, npar, nfiles

    real(dp) :: dim_fract, dim_fract_r,pearson,pearson_r,chi_avr,c, mu, nu

    real(dp), dimension(:,:), allocatable :: X

    real(dp) :: d(Num_dat)
    integer :: N(Num_dat)

    integer, allocatable :: kpar(:)
    real(dp), allocatable :: mupar(:), nupar(:)

    integer :: i, k, kk, ios
    logical :: existe

    character(len=256) :: fichero, nombre, linea

    dim = 2

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
     
        call calc_d(X,d_min,d_max,Num_dat,d,N,P,dim,dim_fract,c,chi,pearson)
        call Calc_d_rad(X,d,Num_dat,P,Pr,d_min,d_max,dim,dim_fract_r,c,pearson_r,chi_avr)

        ! N frente a d de esta trayectoria
        write(nombre,'(A,I0,A)') "dimension/dimension_", k, ".dat"
        open(unit=1, file=trim(nombre), status="replace")
        write(1,'(A)') "# d   N"
        do i = 1, Num_dat
            write(1,*) d(i), N(i)
        end do
        close(1)

        write(2,*) k, mu, nu, dim_fract,dim_fract_r
        write(*,*) "Trayectoria", k, " dimension fractal:", dim_fract, " dimension fractal:", dim_fract_r

        nfiles = nfiles + 1

        deallocate(X)
    end do

    close(2)

    if (nfiles == 0) then
        write(*,*) "No se ha encontrado ninguna trayectoria"
    else
        write(*,*) "Trayectorias analizadas:", nfiles
    end if

    contains
    !Subrutinas de calculo de dimension fractal con cajas
    subroutine Calc_d_Cajas(X,d_min,d_max,Num_dat,d,N,P,dim,dim_fract,c,chi,pearson)
        Implicit none
        integer, intent(in) :: P, dim
        real(dp), intent(in)  :: X(P,dim)
        real(dp), intent(in) :: d_min,d_max
        integer, intent(in) :: Num_dat

        real(dp), intent(out) :: d(Num_dat)
        integer, intent(out) :: N(Num_dat)
        real(dp),intent(out) :: dim_fract,c,chi,pearson
        real(dp) :: sumyy, y_fit
        integer :: n_valid

        real(dp) :: ld(Num_dat),lN(Num_dat), sumx,sumy,sumxy,sumxx
        integer :: i

        call Calc_N(X,d_min,d_max,Num_dat,d,N,P,dim)

        sumx=0.0_dp
        sumy=0.0_dp
        sumxx=0.0_dp
        sumxy=0.0_dp
        n_valid = 0

        do i = 1, Num_dat
            if (N(i) /= 0) then
                n_valid = n_valid + 1
                lN(i) = log(1.0_dp * N(i))
                ld(i) = log(1.0_dp / d(i))

                sumy = sumy + lN(i)
                sumx = sumx + ld(i)
                sumxx = sumxx + ld(i)**2
                sumyy = sumyy + lN(i)**2
                sumxy = sumxy + ld(i) * lN(i)
            else
                print*, "Algo  ha salido mal porque ha salido N=0"
            end if
        end do

        dim_fract = (n_valid * sumxy - sumx * sumy) / (n_valid * sumxx - sumx * sumx)
        c = (sumy - dim_fract * sumx) / real(n_valid, dp)

        !Pearson
        pearson = (n_valid * sumxy - sumx * sumy) / &
                  sqrt((n_valid * sumxx - sumx**2) * (n_valid * sumyy - sumy**2))

        !chi^2 reducido (como sigma=0 es sin los errores)
        chi = 0.0_dp
        do i = 1, Num_dat
            if (N(i) /= 0) then
                y_fit = dim_fract * ld(i) + c
                ! Si tienes errores asociados a lN(i), divide entre sigma_lN(i)**2
                chi = chi + (lN(i) - y_fit)**2
            end if
        end do
        if (n_valid > 2) then
            chi = chi / real(n_valid - 2, dp)
        else
            chi = 0.0_dp
            print*, "Aviso: no hay suficientes puntos para calcular chi^2 reducido"
        end if

    end subroutine
    
    subroutine Calc_N_cajas(X,d_min,d_max,Num_dat,d,N,P,dim)
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

    subroutine calc_Num_cajas(X, d, N, P, dim)
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

    !Calcular la dimension fractal pero con radios de distinto tramaño
    subroutine Calc_d_rad(X,r,Num_dat,P,Pr,f_lo,f_hi,dim,dim_fract,pearson,chi_avr)
        Implicit none
        !P numero de puntos, dim dimension
        integer, intent(in) :: P, Pr, dim
        !Los puntos en cuestion
        real(dp), intent(in)  :: X(P,dim),f_lo,f_hi
        !numero de datos apra el ajuste
        integer, intent(in) :: Num_dat
   
        ! dimension fractal y 
        real(dp),intent(out) :: dim_fract
        !radios, se guardan los ultimos
        real(dp), intent(out) :: r(P),pearson,chi_avr
        !radios donde se empieza
        real(dp)::r0(dim),chi,pearson_aux,y_fit,r_lo,r_hi

        real(dp) :: ld(Num_dat),lN(Num_dat), sumx,sumy,sumyy,sumxy,sumxx,d_aux,dim_fract_aux,c, ci_avr

        integer ::n_valid, i,j,m,k,jmin, jmax

        pearson=0
        chi_avr=0
        dim_fract=0

        !se calculan las distancias al centro
        do k=1,Pr
            do m=1,dim
                r0(m)=X(max(1, (k*(P/2))/Pr),m)
            end do

            do i = 1, P
                r(i)=0
                do j=1,dim
                    r(i)=r(i)+(X(i,j)-r0(j))**2
                end do 
                r(i)=sqrt(r(i))
            end do
        
            call quick_sort(r, P)
            sumy=0
            sumx=0
            sumyy=0
            sumxx=0
            sumxy=0
            do i=1,Num_dat
                !Asumo que (r(i) approx prop sqrt(t)) 
                r_lo = f_lo * r(P)      ! r(P) = distancia máxima a este centro
                r_hi = f_hi * r(P)

                jmin = 2                 ! r(1)=0 (el propio centro), así que empezamos en 2
                do while (jmin < P)
                    if (r(jmin) >= r_lo) exit
                    jmin = jmin + 1
                end do

                jmax = P
                do while (jmax > jmin)
                    if (r(jmax) <= r_hi) exit
                    jmax = jmax - 1
                end do

                j = nint( jmin * (real(jmax,dp)/real(jmin,dp))**( real(i-1,dp)/real(Num_dat-1,dp) ) )
                j = min(max(j, 2), P)
                
                lN(i)=log(1.0_dp*j)
                ld(i)=log(1.0_dp*r(j))

                sumy=sumy+lN(i)
                sumx=sumx+ld(i)
                sumxx=sumxx+ld(i)**2
                sumyy = sumyy + lN(i)**2
                sumxy=sumxy+ld(i)*lN(i)
            end do
             n_valid = Num_dat

            dim_fract_aux=(Num_dat*sumxy-sumx*sumy)*1.0_dp/(Num_dat*sumxx-sumx*sumx)
            dim_fract=dim_fract+dim_fract_aux
            c = (sumy - dim_fract_aux*sumx)/Num_dat

            !Pearson
            pearson_aux = (n_valid * sumxy - sumx * sumy) / sqrt((n_valid * sumxx - sumx**2) * (n_valid * sumyy - sumy**2))
            pearson=pearson+pearson_aux

            !chi^2 reducido (como sigma=0 es sin los errores)
            chi = 0.0_dp
            do i = 1, Num_dat
                y_fit = dim_fract_aux * ld(i) + c
                chi = chi + (lN(i) - y_fit)**2
            end do
            if (n_valid > 2) then
                chi = chi / real(n_valid - 2, dp)
            else
                chi = 0.0_dp
                print*, "Aviso: no hay suficientes puntos para chi^2"
            end if
            chi_avr=chi_avr+chi

        end do
        pearson=pearson/(1.0_dp*Pr)
        chi_avr=chi_avr/(1.0_dp*Pr)
        dim_fract=dim_fract/(1.0_dp*Pr)
    end subroutine

end program