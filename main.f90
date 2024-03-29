
    program main

    use dop853_module
    use dop853_constants

    use easydop853_module

    implicit none

    integer,parameter :: len_data_1961 = 152
    real(wp),dimension(len_data_1961) :: t_1961
    real(wp),dimension(3,len_data_1961) :: pxyz_1961, vxyz_1961

    integer,parameter :: n_max = 100
    real(wp),dimension(1:n_max,0:n_max-1) :: a
    real(wp),dimension(1:n_max) :: b
    real(wp),dimension(2:n_max,0:n_max-2) :: c
    real(wp),dimension(2:n_max) :: d
    real(wp),dimension(2:n_max,0:n_max-1) :: e
    real(wp),dimension(0:n_max,0:n_max) :: p_matrix
    real(wp),dimension(0:n_max,0:n_max) :: c_matrix, s_matrix

    real(wp),parameter :: s = sqrt(6378136.3D0**3/3986004.415D8)

    integer               :: i
    real(wp)              :: t
    real(wp),dimension(6) :: pv
    real(wp)              :: a_echo, e_echo

    call readecho(len_data_1961, t_1961, pxyz_1961, vxyz_1961)

    call getabcdep(n_max, a, b, c, d, e, p_matrix)
    call readegm(n_max, c_matrix, s_matrix)

    pv(1:3) = pxyz_1961(:,1)
    pv(4:6) = vxyz_1961(:,1)
    !call pv2ae(pv(1:3), pv(4:6), a_echo, e_echo)
    !call pv2m(pv(1:3), pv(4:6))
    !print *, a_echo*(1-e_echo)*6.3781363_wp, e_echo
    do i = 1, len_data_1961-1
        t = (i-1)*86400/s
        !pv(1:3) = pxyz_1961(:,i+1)
        !pv(4:6) = vxyz_1961(:,i+1)
        !call easydop853(f_nosrp, t, i*86400/s, pv)
        !call easydop853(f_srp1, t, i*86400/s, pv)
        !call easydop853(f_srp2, t, i*86400/s, pv)
        call easydop853(f_srp3, t, i*86400/s, pv)
        !print *, pv
        call pv2ae(pv(1:3), pv(4:6), a_echo, e_echo)
        !call pv2m(pv(1:3), pv(4:6))
        print *, a_echo*(1-e_echo)*6.3781363_wp, e_echo
        !print *, i+1
        !print *, pxyz_1961(:,i+1)
        !print *, vxyz_1961(:,i+1)
        !call pv2ae(pxyz_1961(:,i+1), vxyz_1961(:,i+1), a_echo, e_echo)
        !call pv2m(pxyz_1961(:,i+1), vxyz_1961(:,i+1))
        !print *, a_echo*(1-e_echo)*6.3781363_wp, e_echo
        !print *, pv(1:3)
        !print *, pv(4:6)
        !call pv2ae(pv(1:3), pv(4:6), a_echo, e_echo)
        !call pv2m(pv(1:3), pv(4:6))
        !print *, a_echo*(1-e_echo)*6.3781363_wp, e_echo
        !print *, ''
    end do

    contains

    subroutine f_nosrp(me,t,pv,vf)

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: pv
    real(wp),dimension(:),intent(out) :: vf

    real(wp),dimension(3) :: fe
    real(wp),dimension(3) :: fp

    call getfe(t_1961(1)+t*s/86400, pv(1:3), &
               n_max, a, b, c, d, e, p_matrix, c_matrix, s_matrix, fe)
    call getfp(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fp)

    vf(1:3) = pv(4:6)
    vf(4:6) = fe+fp
    !vf(4:6) =-pv(1:3)/sum(pv(1:3)**2.0_wp)**1.5_wp

    end subroutine f_nosrp

    subroutine f_srp1(me,t,pv,vf)

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: pv
    real(wp),dimension(:),intent(out) :: vf

    real(wp),dimension(3) :: fe
    real(wp),dimension(3) :: fp
    real(wp),dimension(3) :: fsrp

    call getfe(t_1961(1)+t*s/86400, pv(1:3), &
               n_max, a, b, c, d, e, p_matrix, c_matrix, s_matrix, fe)
    call getfp(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fp)
    call srp1(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fsrp)

    vf(1:3) = pv(4:6)
    vf(4:6) = fe+fp+fsrp

    end subroutine f_srp1

    subroutine f_srp2(me,t,pv,vf)

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: pv
    real(wp),dimension(:),intent(out) :: vf

    real(wp),dimension(3) :: fe
    real(wp),dimension(3) :: fp
    real(wp),dimension(3) :: fsrp

    call getfe(t_1961(1)+t*s/86400, pv(1:3), &
               n_max, a, b, c, d, e, p_matrix, c_matrix, s_matrix, fe)
    call getfp(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fp)
    call srp2(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fsrp)

    vf(1:3) = pv(4:6)
    vf(4:6) = fe+fp+fsrp

    end subroutine f_srp2

    subroutine f_srp3(me,t,pv,vf)

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: pv
    real(wp),dimension(:),intent(out) :: vf

    real(wp),dimension(3) :: fe
    real(wp),dimension(3) :: fp
    real(wp),dimension(3) :: fsrp

    call getfe(t_1961(1)+t*s/86400, pv(1:3), &
               n_max, a, b, c, d, e, p_matrix, c_matrix, s_matrix, fe)
    call getfp(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fp)
    call srp3(2400000.5D0+t_1961(1)+t*s/86400, pv(1:3), fsrp)

    vf(1:3) = pv(4:6)
    vf(4:6) = fe+fp+fsrp

    end subroutine f_srp3

    end program main
