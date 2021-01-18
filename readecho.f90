
    subroutine readecho(len_data_1961,t_1961,xyz_1961)

        use iso_fortran_env,    only: wp => real64

        implicit none

        integer,intent(in) :: len_data_1961
        real(wp),dimension(len_data_1961),intent(out) :: t_1961
        real(wp),dimension(3,len_data_1961),intent(out) :: xyz_1961

        integer :: i
        real(wp),dimension(7,len_data_1961) :: data_1961
        real(wp),dimension(len_data_1961)   :: e, M, q, Ra, ic, ap
        real(wp),dimension(len_data_1961)   :: Ea, a
        real(wp),dimension(len_data_1961)   :: x, y

        open(10,file='data_1961')

        do i = 1, len_data_1961
            read(10,'(f8.2,f7.3,f10.5,f8.5,f8.6,f9.6,f9.7)') data_1961(:,i)
        end do
        print *, data_1961(1,:)
        print *, data_1961(2,:)
        print *, data_1961(3,:)
        print *, data_1961(4,:)
        print *, data_1961(5,:)
        print *, data_1961(6,:)
        print *, data_1961(7,:)
        close(10)

        data_1961(3,:) = data_1961(3,:)+3.508e-5_wp*(data_1961(1,:)-33281.0_wp)

        data_1961(2:4,:) = data_1961(2:4,:)/180.0_wp*3.141592653589793_wp
        data_1961(6,:) = data_1961(6,:)*3.141592653589793_wp
        data_1961(7,:) = data_1961(7,:)*1.0e6_wp

        t_1961 = data_1961(1,:)

        e = data_1961(5,:)
        M = data_1961(6,:)
        q = data_1961(7,:)

        Ea = M
        if (sum((M+e*sin(Ea)-Ea)**2.0_wp).ge.1.0e-16_wp) then
            Ea = M+e*sin(Ea)
        end if

        a = q/(1.0_wp-e)
        x = a*(cos(Ea)-e)
        y = a*(sqrt(1.0_wp-e**2.0_wp)*sin(Ea))
        
        Ra = data_1961(3,:)
        ic = data_1961(4,:)
        ap = data_1961(2,:)

        xyz_1961(1,:) = (cos(Ra)*cos(ap)-sin(Ra)*cos(ic)*sin(ap))*x &
                       -(cos(Ra)*sin(ap)+sin(Ra)*cos(ic)*cos(ap))*y
        xyz_1961(2,:) = (sin(Ra)*cos(ap)+cos(Ra)*cos(ic)*sin(ap))*x &
                       -(sin(Ra)*sin(ap)-cos(Ra)*cos(ic)*cos(ap))*y
        xyz_1961(3,:) = sin(ic)*sin(ap)*x+sin(ic)*cos(ap)*y

    end subroutine readecho
