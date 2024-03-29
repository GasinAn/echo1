
    subroutine readecho(len_data_1961,t_1961,pxyz_1961,vxyz_1961)

        use iso_fortran_env,    only: wp => real64

        implicit none

        integer,intent(in) :: len_data_1961
        real(wp),dimension(len_data_1961),intent(out) :: t_1961
        real(wp),dimension(3,len_data_1961),intent(out) :: pxyz_1961, vxyz_1961

        integer                             :: i
        double precision                    :: t2m
        double precision                    :: RMATN(3,3)
        double precision                    :: eo
        double precision                    :: iau_EO06A
        real(wp)                            :: tt, ut
        real(wp),dimension(7,len_data_1961) :: data_1961
        real(wp),dimension(len_data_1961)   :: e, M, q, Ra, ic, ap
        real(wp),dimension(len_data_1961)   :: Ea, a
        real(wp),dimension(len_data_1961)   :: px, py, vx, vy
        real(wp),dimension(3,3)             :: rc2i, ri2c

        open(10,file='data_1961')

        do i = 1, len_data_1961
            read(10,'(f8.2,f7.3,f10.5,f8.5,f8.6,f9.6,f9.7)') data_1961(:,i)
        end do

        close(10)

        data_1961(3,:) = data_1961(3,:)+3.508e-5_wp*(data_1961(1,:)-33281.0_wp)

        data_1961(2:4,:) = data_1961(2:4,:)/180.0_wp*3.141592653589793_wp
        data_1961(6,:) = data_1961(6,:)*2.0_wp*3.141592653589793_wp
        data_1961(7,:) = data_1961(7,:)*1.0e6_wp

        t_1961 = data_1961(1,:)

        e = data_1961(5,:)
        M = data_1961(6,:)
        q = data_1961(7,:)

        Ea = M
        if (sum((M+e*sin(Ea)-Ea)**2.0_wp).ge.1.0e-16_wp) then
            Ea = M+e*sin(Ea)
        end if

        a = q/(1.0_wp-e)/6378136.3_wp
        px = a*(+cos(Ea)-e)
        py = a*(sqrt(1.0_wp-e**2.0_wp)*sin(Ea))
        vx = sqrt(a/(px**2.0_wp+py**2.0_wp))*(-sin(Ea))
        vy = sqrt(a/(px**2.0_wp+py**2.0_wp))*(sqrt(1.0_wp-e**2.0_wp)*cos(Ea))

        Ra = data_1961(3,:)
        ic = data_1961(4,:)
        ap = data_1961(2,:)

        pxyz_1961(1,:) = (cos(Ra)*cos(ap)-sin(Ra)*cos(ic)*sin(ap))*px &
                        -(cos(Ra)*sin(ap)+sin(Ra)*cos(ic)*cos(ap))*py
        pxyz_1961(2,:) = (sin(Ra)*cos(ap)+cos(Ra)*cos(ic)*sin(ap))*px &
                        -(sin(Ra)*sin(ap)-cos(Ra)*cos(ic)*cos(ap))*py
        pxyz_1961(3,:) =  sin(ic)*sin(ap)*px+sin(ic)*cos(ap)*py

        vxyz_1961(1,:) = (cos(Ra)*cos(ap)-sin(Ra)*cos(ic)*sin(ap))*vx &
                        -(cos(Ra)*sin(ap)+sin(Ra)*cos(ic)*cos(ap))*vy
        vxyz_1961(2,:) = (sin(Ra)*cos(ap)+cos(Ra)*cos(ic)*sin(ap))*vx &
                        -(sin(Ra)*sin(ap)-cos(Ra)*cos(ic)*cos(ap))*vy
        vxyz_1961(3,:) =  sin(ic)*sin(ap)*vx+sin(ic)*cos(ap)*vy

        ! V(true) = RMATN * V(mean)
        ! [?    [... ... ...    [1
        !  ?  =  ... ... ...  *  0
        !  ?]    ... ... ...]    0]
        ! [?    [RMATN(1, 1)
        !  ?  =  RMATN(2, 1)
        !  ?]    RMATN(3, 1)]
        ! true -> mean = atan2(RMATN(2, 1), RMATN(1, 1))
        ! mean -> true = -atan2(RMATN(2, 1), RMATN(1, 1))
        do i = 1, len_data_1961
            tt = t_1961(i)
            call iau_NUM06A(2400000.5D0, tt, RMATN)
            t2m = atan2(RMATN(2, 1), RMATN(1, 1))
            pxyz_1961(1:2,i) = [cos(-t2m)*pxyz_1961(1,i) &
                               +sin(-t2m)*pxyz_1961(2,i), &
                               -sin(-t2m)*pxyz_1961(1,i) &
                               +cos(-t2m)*pxyz_1961(2,i)]
            vxyz_1961(1:2,i) = [cos(-t2m)*vxyz_1961(1,i) &
                               +sin(-t2m)*vxyz_1961(2,i), &
                               -sin(-t2m)*vxyz_1961(1,i) &
                               +cos(-t2m)*vxyz_1961(2,i)]
        end do

        ! cio -> eqx = eo
        ! eqx -> cio = -eo
        do i = 1, len_data_1961
            tt = t_1961(i)
            eo = iau_EO06A(2400000.5D0, tt)
            pxyz_1961(1:2,i) = [cos(-eo)*pxyz_1961(1,i) &
                               +sin(-eo)*pxyz_1961(2,i), &
                               -sin(-eo)*pxyz_1961(1,i) &
                               +cos(-eo)*pxyz_1961(2,i)]
            vxyz_1961(1:2,i) = [cos(-eo)*vxyz_1961(1,i) &
                               +sin(-eo)*vxyz_1961(2,i), &
                               -sin(-eo)*vxyz_1961(1,i) &
                               +cos(-eo)*vxyz_1961(2,i)]
        end do

        do i = 1, len_data_1961
            call iau_C2I06A(2400000.5D0, t_1961(i), rc2i)
            ri2c = reshape( &
               [rc2i(2,2)*rc2i(3,3)-rc2i(2,3)*rc2i(3,2), &
                rc2i(2,3)*rc2i(3,1)-rc2i(2,1)*rc2i(3,3), &
                rc2i(2,1)*rc2i(3,2)-rc2i(2,2)*rc2i(3,1), &
                rc2i(3,2)*rc2i(1,3)-rc2i(3,3)*rc2i(1,2), &
                rc2i(3,3)*rc2i(1,1)-rc2i(3,1)*rc2i(1,3), &
                rc2i(3,1)*rc2i(1,2)-rc2i(3,2)*rc2i(1,1), &
                rc2i(1,2)*rc2i(2,3)-rc2i(1,3)*rc2i(2,2), &
                rc2i(1,3)*rc2i(2,1)-rc2i(1,1)*rc2i(2,3), &
                rc2i(1,1)*rc2i(2,2)-rc2i(1,2)*rc2i(2,1)], &
                [3,3])
            pxyz_1961(1:3,i) = [ri2c(1,1)*pxyz_1961(1,i) &
                               +ri2c(1,2)*pxyz_1961(2,i) &
                               +ri2c(1,3)*pxyz_1961(3,i), &
                                ri2c(2,1)*pxyz_1961(1,i) &
                               +ri2c(2,2)*pxyz_1961(2,i) &
                               +ri2c(2,3)*pxyz_1961(3,i), &
                                ri2c(3,1)*pxyz_1961(1,i) &
                               +ri2c(3,2)*pxyz_1961(2,i) &
                               +ri2c(3,3)*pxyz_1961(3,i)]
            vxyz_1961(1:3,i) = [ri2c(1,1)*vxyz_1961(1,i) &
                               +ri2c(1,2)*vxyz_1961(2,i) &
                               +ri2c(1,3)*vxyz_1961(3,i), &
                                ri2c(2,1)*vxyz_1961(1,i) &
                               +ri2c(2,2)*vxyz_1961(2,i) &
                               +ri2c(2,3)*vxyz_1961(3,i), &
                                ri2c(3,1)*vxyz_1961(1,i) &
                               +ri2c(3,2)*vxyz_1961(2,i) &
                               +ri2c(3,3)*vxyz_1961(3,i)]
        end do

    end subroutine readecho
