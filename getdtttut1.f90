
    subroutine getdtttut1(ttmjd, dtttut1mjd)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: dtttut1mjd

        dtttut1mjd = -8.44128673D-13*ttmjd**2 &
                     +1.01483271D-07*ttmjd    &
                     -2.25794010D-03          &
                     +7.91984037D-06*sin(8.95765509D-04*ttmjd-5.18352618D-01)

    end subroutine getdtttut1
