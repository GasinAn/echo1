
    subroutine getpmy(ttmjd, pmy)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: pmy

        pmy = -1.28468440D-01*sin(1.45140969D-02*ttmjd+4.74212791D+00) &
              -8.69726805D-02*sin(1.72063784D-02*ttmjd+1.02303508D+00) &
              +8.48529906D-06*ttmjd-1.13657698D-01

    end subroutine getpmy
