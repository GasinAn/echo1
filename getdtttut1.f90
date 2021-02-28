
    subroutine getdtttut1(ttmjd, dtttut1mjd)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: dtttut1mjd

        dtttut1mjd = -8.44142565D-13*ttmjd**2 &
                     +1.01484504D-07*ttmjd    &
                     -2.25796675D-03          &
                     +7.91932922D-06*sin(8.95775798D-04*ttmjd-5.18784739D-01)

    end subroutine getdtttut1
