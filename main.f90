
    program main
!-------------------------------------------------------------------------------
    use dop853_module
    use dop853_constants

    use easydop853_module
!-------------------------------------------------------------------------------
    implicit none

    integer,parameter :: n_max = 100
    real(wp),dimension(0:n_max,0:n_max) :: c, s

    integer,parameter :: len_data_1961 = 110
    real(wp),dimension(len_data_1961) :: t_1961
    real(wp),dimension(3,len_data_1961) :: xyz_1961

    procedure(deriv_func) :: fcn
    !! subroutine computing the value of \(dy/dx=f(x,y)\)
    real(wp)              :: x  = 0.0_wp
    !! `x` value (input is initial value and output is final value)
    real(wp)              :: xf = 225000.0_wp
    !! endpoint of integration (final value of `x`)
    real(wp),dimension(6) :: y  = [1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp]
    !! `y` value (input is initial value and output is final value)    
!-------------------------------------------------------------------------------
    call readegm(n_max,c,s)

    call readecho(len_data_1961,t_1961,xyz_1961)

    call easydop853(fcn,x,xf,y)
    
    print *, x
    print *, y

    end program main
