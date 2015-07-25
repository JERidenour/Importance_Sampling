program impSamp
        interface
                function randU(d)
                        integer, intent(in) :: d
                        real, dimension(d) :: randU
                        integer, parameter :: d15 = selected_real_kind(15)
                        integer :: time_int
                        real(kind = d15) :: x, m, time, a
                        integer, dimension(8) :: values
                end function randU
                function f_func(x,d)
                        integer, intent(in) :: d
                        real, dimension(d), intent(in) :: x
                        real, dimension(d) :: f_func
                        integer, parameter :: d15 = selected_real_kind(15)
                        real(kind = d15) :: pi
                end function f_func
                function g_func(x,d)
                        integer, intent(in) :: d
                        real, dimension(d), intent(in) :: x
                        real, dimension(d) :: g_func
                        integer, parameter :: d15 = selected_real_kind(15)
                        real(kind = d15) :: pi
                end function g_func
                function phi_func(x,d)
                        integer, intent(in) :: d
                        real, dimension(d), intent(in) :: x
                        real, dimension(d) :: phi_func
                end function phi_func
                function w_func(x,d)
                        integer, intent(in) :: d
                        real, dimension(d), intent(in) :: x
                        real, dimension(d) :: w_func, f, g
                end function w_func
                function stdn(d)
                        integer, intent(in) :: d
                        real, dimension(d) :: stdn
                        real :: z1, z2, u1, u2
                        integer, parameter :: d15 = selected_real_kind(15)
                        real(kind = d15) :: pi
                end function stdn
        end interface 

        integer :: d=50000
        real, dimension(50000) :: X, phi, w, prod
        real :: tau
        print *, "Importance Sampling for estimating the following integral.."
        print *, "integrate from -0.5 to 0.25 exp(sin(x))*1/(pi*(1 + x*x)) dx"

        !draw from g (the standard normal)
        X = stdn(d)
        
        !create importance weights 
        w = w_func(X,d)

        !estimate phi*w with standard MC
        phi = phi_func(X,d)
        do i=1,d
                prod(i) = phi(i)*w(i)/d
        end do
        tau = sum(prod)

        !the expectation value is now an approximation of the integral
        print *, "The integral value (50000 samples):"
        print *, tau
end program impSamp

!multiplicative congruential generator, creates a vector
!of size d with uniformly distributed psuedorandom numbers in (0,1)
function randU(d)
        integer, intent(in) :: d
        real, dimension(d) :: randU
        integer, parameter :: d15 = selected_real_kind(15)
        integer :: time_int
        real(kind = d15) :: x, m, time, a
        integer, dimension(8) :: values
        call date_and_time(VALUES=values)
        m = 4294967295_d15
        time_int = values(8)
        time = real(time_int)
        a = 16807
        x = mod(a*time,m)
        do i=1,d
                randU(i) = x/m
                x = mod(a*x,m)
        end do
end function randU

!the function you want to integrate
!returns the function value of exp(sin(x)) * 1/(pi*(1+x^2)) for 
!input values x
function f_func(x,d)
        integer, intent(in) :: d
        real, dimension(d), intent(in) :: x
        real, dimension(d) :: f_func
        integer, parameter :: d15 = selected_real_kind(15)
        real(kind = d15) :: pi
        pi = 3.14159265359_d15
        do i=1,d
                f_func(i) = exp(sin(x(i)))*1/(pi*(1 + x(i)*x(i)))
        end do
end function f_func

!the density function of the intrumental distribution (here, the 
!standard normal)
!returnes the function value of 1/sqrt(2*pi)*exp(-0.5*x^2) for
!input values x 
function g_func(x,d)
        integer, intent(in) :: d
        real, dimension(d), intent(in) :: x
        real, dimension(d) :: g_func
        integer, parameter :: d15 = selected_real_kind(15)
        real(kind = d15) :: pi
        pi = 3.14159265359_d15
        do i=1,d
                g_func(i) = 1/sqrt(2*pi)*exp(-0.5*x(i)*x(i))
        end do
end function g_func

!returnes an array with elements 1 if the corresponding value of 
!x is in [-0.5 0.25] and 0 otherwise
function phi_func(x,d)
        integer, intent(in) :: d
        real, dimension(d), intent(in) :: x
        real, dimension(d) :: phi_func
        do i=1,d
                if(x(i) >= -0.5 .and. x(i) <= 0.25) then
                        phi_func(i) = 1
                else
                        phi_func(i) = 0
                end if
        end do
end function phi_func

!the weight function 
!defined as f/g
function w_func(x,d)
        interface
                function f_func(x,d)
                        integer, intent(in) :: d
                        real, dimension(d), intent(in) :: x
                        real, dimension(d) :: f_func
                        integer, parameter :: d15 = selected_real_kind(15)
                        real(kind = d15) :: pi
                end function f_func
                function g_func(x,d)
                        integer, intent(in) :: d
                        real, dimension(d), intent(in) :: x
                        real, dimension(d) :: g_func
                        integer, parameter :: d15 = selected_real_kind(15)
                        real(kind = d15) :: pi
                end function g_func
        end interface
        integer, intent(in) :: d
        real, dimension(d), intent(in) :: x
        real, dimension(d) :: w_func, f, g
        f = f_func(x,d)
        g = g_func(x,d)
        do i=1,d
                w_func(i) = f(i)/g(i)
        end do
end function w_func

!creates d draws from a standard normal distribution
!uses the Box–Muller transform 
function stdn(d)
        interface
                function randU(d)
                        integer, intent(in) :: d
                        real, dimension(d) :: randU
                        integer, parameter :: d15 = selected_real_kind(15)
                        integer :: time_int
                        real(kind = d15) :: x, m, time, a
                        integer, dimension(8) :: values
                end function randU
        end interface
        integer, intent(in) :: d
        real, dimension(d) :: stdn, u
        real :: z1, z2
        integer, parameter :: d15 = selected_real_kind(15)
        real(kind = d15) :: pi

        pi = 3.14159265359_d15
        u = randU(d+1)
        i = 1
        do while (i<=d)
                z1 = sqrt(-2*log(u(i)))*cos(2*pi*u(i+1))
                stdn(i) = z1
                i = i+1
                z2 = sqrt(-2*log(u(i)))*sin(2*pi*u(i+1))
                stdn(i) = z2
        end do
end function
