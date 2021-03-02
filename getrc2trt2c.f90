
    subroutine getrc2trt2c(ttmjd, rc2t, rt2c)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),dimension(3,3),intent(out) :: rc2t
        real(wp),dimension(3,3),intent(out) :: rt2c

        real(wp) :: utmjd
        real(wp) :: xp
        real(wp) :: yp

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

    end subroutine getrc2trt2c

