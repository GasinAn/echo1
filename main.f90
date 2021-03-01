
    program main

    use dop853_module
    use dop853_constants

    use easydop853_module

    implicit none

    integer,parameter :: len_data_1961 = 152
    real(wp),dimension(len_data_1961) :: t_1961
    real(wp),dimension(3,len_data_1961) :: pxyz_1961, vxyz_1961

    real(wp),parameter :: s = sqrt(6378136.3D0**3/3986004.415D8)

    procedure(deriv_func) :: fcn
    real(wp)              :: t
    real(wp),dimension(6) :: pv    

    call readecho(len_data_1961, t_1961, pxyz_1961, vxyz_1961)

    t = t_1961(1)
    pv(1:3) = pxyz_1961(:,1)
    pv(4:6) = vxyz_1961(:,1)
    print *, pv
    call easydop853(fcn, t, t+(len_data_1961-1)*86400/s, pv)
    print *, pv

    print *, pxyz_1961(:,len_data_1961)
    print *, vxyz_1961(:,len_data_1961)

    end program main
