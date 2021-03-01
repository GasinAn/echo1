
    program main
!-------------------------------------------------------------------------------
    use dop853_module
    use dop853_constants

    use easydop853_module
!-------------------------------------------------------------------------------
    implicit none

    integer,parameter :: len_data_1961 = 152
    real(wp),dimension(len_data_1961) :: t_1961
    real(wp),dimension(3,len_data_1961) :: pxyz_1961, vxyz_1961

    integer,parameter :: n_max = 100
    real(wp),dimension(0:n_max,0:n_max) :: c_matrix, s_matrix

    real(wp),dimension(1:n_max,0:n_max-1) :: a
    real(wp),dimension(1:n_max) :: b
    real(wp),dimension(2:n_max,0:n_max-2) :: c
    real(wp),dimension(2:n_max) :: d
    real(wp),dimension(2:n_max,0:n_max-1) :: e
    real(wp),dimension(0:n_max,0:n_max) :: p
    real(wp),dimension(3) :: f

    procedure(deriv_func) :: fcn
    !! subroutine computing the value of \(dy/dx=f(x,y)\)
    real(wp)              :: x  = 0.0_wp
    !! `x` value (input is initial value and output is final value)
    real(wp)              :: xf = 225000.0_wp
    !! endpoint of integration (final value of `x`)
    real(wp),dimension(6) :: y  = [1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp]
    !! `y` value (input is initial value and output is final value)    
!-------------------------------------------------------------------------------
    call readecho(len_data_1961, t_1961, pxyz_1961, vxyz_1961)

    print *, 'readecho ok!'

    call readegm(n_max, c_matrix, s_matrix)

    print *, 'readegm ok!'

    call getabcdep(n_max, a, b, c, d, e, p)

    print *, 'getabcdep ok!'

    call getfe(36294.0D0, 2.0_wp**(1.0_wp/3.0_wp), 0.0_wp, 0.0_wp, &
               n_max, a, b, c, d, e, c_matrix, s_matrix, p, f)

    print *, sum(f**2.0_wp)**0.5_wp

    call getfp(2436294.5_wp, [8000000.0_wp/6378136.3_wp,0.0_wp,0.0_wp], f)

    print *, sum(f**2.0_wp)**0.5_wp

    call easydop853(fcn, x, xf, y)

    print *, x
    print *, y

    end program main
