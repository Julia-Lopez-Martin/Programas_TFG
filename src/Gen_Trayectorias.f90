Program Generar_trayectorias
    !gfortran -O3 Generar_trayectorias.f90 -o Generar_trayectorias.exe && Generar_trayectorias.exe
    use rutas
    implicit none
    ! ---------- SELECCION DE MODELO ----------
    integer, parameter :: modelo = 1

    ! ---------- PARAMETROS GENERALES ----------
    integer, parameter :: dp = selected_real_kind(15,300)
    integer, parameter :: N=10000000
    integer, parameter :: dim=3
    real(dp), parameter :: PI=3.141592653589793238462_dp

    ! ---------- PARAMETROS A MODIFICAR (modelo 1) ----------
    integer, parameter :: n_med = 5
    real(dp), parameter :: mu_min = 10.0_dp
    real(dp), parameter :: mu_max = 10.0_dp
    real(dp), parameter :: nu_min = 10.0_dp
    real(dp), parameter :: nu_max = 10.0_dp

    integer :: i,j,k,b_m
    !parametros auxiliares para langevine
    real(dp) :: v_aux(dim),x_aux(dim),mu,nu,dt,d

    character(len=100) :: nombre
    integer :: ios
    call random_seed()

    dt=0.01_dp
    b_m=100
    
    !Hace el directorio de trayectorias donde se guardaran los datos generados
    call execute_command_line("mkdir -p " // DIR_TRAY)

    !hace el archivo de aprametros con todos los parametros comprueba si se ha creado correctamente y si no da error
    open(newunit=u_par, file=DIR_TRAY // "parametros.dat", status="replace", iostat=ios)
    if (ios /= 0) stop "No se pudo crear parametros.dat"

    !escribe los parametros
   
    write(u_par,*) "# Parametros de la ejecucion"
    write(u_par,*) "# N =", N, " b_m =", b_m, " modelo =", modelo
    write(u_par,*) "# num mu nu"

    do k=1,n_med
        ! Valores de mu y nu, ajustar si no se quieren equiespaciados
        mu = mu_min + (k-1)*(mu_max-mu_min)/real(n_med-1,dp)
        nu = nu_min + (k-1)*(nu_max-nu_min)/real(n_med-1,dp)

        !Inicializar todas las variables
        do i=1,dim
            x_aux(i)=0
            v_aux(i)=0
        end do
        d=0

        write(nombre,'(A,I0,A)') DIR_TRAY // "trayectoria_", k, ".dat"
        open(newunit=u_tray, file=trim(nombre), status="replace", access="stream", form="unformatted")
        if (ios /= 0) stop "No se pudo abrir el archivo de trayectoria"

        do i=0,N
            do j=0,b_m
                select case (modelo)
                    case (1)
                    call Langevine_next_step(x_aux,v_aux,dt,mu,nu,d)
                    case (2)
                    ! pendiente
                    case default
                    stop "modelo debe ser 1 o 2"
                end select
            end do

            write(u_tray) x_aux(1:dim)

        end do
        
        select case (modelo)
            case (1)
                write(u_par,*) k, mu, nu
                !Indicación en pantalla de la trayectoria
                write(*,*) "# trayectoria:", k 
                write(*,*) "mu:", mu,"nu:", nu,"x_final:", x_aux(1:dim)
                write(*,*) "sqrtx^2_promedio:", sqrt(d/((real(N+1,dp)*real(b_m+1,dp))))
        end select
    end do
    close(u_par)

    contains
    ! ---------------- MODELO 1 ----------------
    Subroutine Langevine_next_step(x,v,dt,mu,nu,d,dim)
        implicit none
        real(dp), intent(in) :: dt, mu, nu
        integer, intent(in) :: dim
        real(dp), intent(inout) :: x(dim),v(dim),d
        real(dp) :: ruido, u,u2,d_aux
        integer :: i

        do i=1,dim
            call random_number(u)
            call random_number(u2)
            u=1.0_dp - u
            ruido=sqrt(-2.0_dp*log(u))*cos(2*PI*u2)
            v(i)= v(i) - mu*v(i)*dt + nu*sqrt(dt)*ruido
            d_aux=x(i)
            x(i)=x(i)+dt*v(i)
            d_aux=x(i)-d_aux
            d=d+d_aux**2
        end do
    end subroutine

end program