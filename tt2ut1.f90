
    subroutine tt2ut1(ttmjd, ut1mjd)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: ut1mjd

        real(wp) :: y, dy, dt

        y = 1961+(ttmjd-(37300+365.25_wp/24))/365.25_wp

        dy = y-1950

        dt = 29.07_wp+0.407_wp*dy-4.2918e-3_wp*dy**3+3.926187e-4*dy**4
        ut1mjd = ttmjd+dt/86400
        !ut1mjd = ttmjd

    end subroutine tt2ut1
