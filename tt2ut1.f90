
    subroutine tt2ut1(ttmjd, ut1mjd)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: ut1mjd

        real(wp) :: y, dy, dt

        !y = 1961+(ttmjd-37300)/365.25

        !dy = y-1950
        !dy = y-1975

        !dt = 29.07+0.407*dy-dy**2/223+dy**3/2547
        !dt = 45.45+1.067*dy-dy**2/260-dy**3/718
        dt = 33.2_wp+(35.7_wp-33.2_wp)*((ttmjd-36934)/(38761-36934))

        ut1mjd = ttmjd-dt/86400
        !ut1mjd = ttmjd-(dt+0.2_wp)/86400
        !ut1mjd = ttmjd-(dt-0.2_wp)/86400
        !ut1mjd = ttmjd-33.2_wp/86400
        !ut1mjd = ttmjd-35.7_wp/86400
        !ut1mjd = ttmjd

    end subroutine tt2ut1
