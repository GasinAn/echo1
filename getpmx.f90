
    subroutine getpmx(ttmjd, pmx)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: ttmjd
        real(wp),intent(out) :: pmx

        pmx = -1.29058966D-01*sin(1.45170623D-02*ttmjd-3.26630385D+00) &
              +9.62891298D-02*sin(1.72036211D-02*ttmjd+2.75015601D+00) &
              +1.56778978D-10*ttmjd**2-9.47440461D-06*ttmjd+1.28111560D-01
        pmx = pmx/3600/180*3.141592653589793_wp
        !pmx = 0.0_wp

    end subroutine getpmx
