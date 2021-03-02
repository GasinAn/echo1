
    program main

    use dop853_module
    use dop853_constants

    use easydop853_module

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

    real(wp),parameter :: s = sqrt(6378136.3D0**3/3986004.415D8)

    real(wp)              :: t
    real(wp),dimension(6) :: pv    

    call readecho(len_data_1961, t_1961, pxyz_1961, vxyz_1961)

    call readegm(n_max, c_matrix, s_matrix)
    call getabcdep(n_max, a, b, c, d, e, p)

    t = 0
    pv(1:3) = pxyz_1961(:,1)
    pv(4:6) = vxyz_1961(:,1)
    call easydop853(fnosrp, t, 86400/s, pv)
    print *, pv

    print *, pxyz_1961(:,2)
    print *, vxyz_1961(:,2)

    contains

    subroutine fnosrp(me,t,pv,vf)

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: pv
    real(wp),dimension(:),intent(out) :: vf

    real(wp),dimension(3) :: fe
    real(wp),dimension(3) :: fp

    call getfe(t_1961(1)+t*s/86400, pv(1), pv(2), pv(3), &
               n_max, a, b, c, d, e, c_matrix, s_matrix, p, fe)
    call getfp(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fp)

    vf(1:3) = pv(4:6)
    vf(4:6) = fe+fp

    end subroutine fnosrp

    end program main
