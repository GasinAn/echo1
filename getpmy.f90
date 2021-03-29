
    subroutine getpmy(ttmjd, pmy)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: pmy

        pmy = -1.28251299D-01*sin(1.45147961D-02*ttmjd-1.57508005D+00) &
              -8.69414322D-02*sin(1.72055068D-02*ttmjd+1.06600427D+00) &
              -4.07162241D-10*ttmjd**2+4.79325440D-05*ttmjd-1.05335310D+00
        pmy = pmy/3600/180*3.141592653589793_wp
        !pmy = 0.0_wp

    end subroutine getpmy
