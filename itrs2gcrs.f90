
    subroutine itrs2gcrs(ttmjd, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),dimension(3),intent(inout) :: f

        real(wp) :: utmjd
        real(wp) :: xp
        real(wp) :: yp
        real(wp),dimension(3,3) :: rc2t
        real(wp),dimension(3,3) :: rt2c

        call tt2ut1(ttmjd, utmjd)
        call getpmx(ttmjd, xp)
        call getpmy(ttmjd, yp)

        call iau_C2T06A(2400000.5D0, ttmjd, 2400000.5D0, utmjd, xp, yp, rc2t)

        rt2c = reshape( &
               [rc2t(2,2)*rc2t(3,3)-rc2t(2,3)*rc2t(3,2), &
                rc2t(2,3)*rc2t(3,1)-rc2t(2,1)*rc2t(3,3), &
                rc2t(2,1)*rc2t(3,2)-rc2t(2,2)*rc2t(3,1), &
                rc2t(3,2)*rc2t(1,3)-rc2t(3,3)*rc2t(1,2), &
                rc2t(3,3)*rc2t(1,1)-rc2t(3,1)*rc2t(1,3), &
                rc2t(3,1)*rc2t(1,2)-rc2t(3,2)*rc2t(1,1), &
                rc2t(1,2)*rc2t(2,3)-rc2t(1,3)*rc2t(2,2), &
                rc2t(1,3)*rc2t(2,1)-rc2t(1,1)*rc2t(2,3), &
                rc2t(1,1)*rc2t(2,2)-rc2t(1,2)*rc2t(2,1)], &
                [3,3])

        f = [rt2c(1,1)*f(1)+rt2c(1,2)*f(2)+rt2c(1,3)*f(3), &
             rt2c(2,1)*f(1)+rt2c(2,2)*f(2)+rt2c(2,3)*f(3), &
             rt2c(3,1)*f(1)+rt2c(3,2)*f(2)+rt2c(3,3)*f(3)]

    end subroutine itrs2gcrs
