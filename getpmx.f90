
    subroutine getpmx(ttmjd, pmx)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: pmx

        pmx = -1.29111174D-01*sin(1.45167959D-02*ttmjd+3.02874420D+00) &
              +9.62902744D-02*sin(1.72034659D-02*ttmjd+2.75757541D+00) &
              +5.71426777D-06*ttmjd-2.33690344D-01
        pmx = pmx/3600/180*3.141592653589793_wp

    end subroutine getpmx
