
    subroutine srp1(TDB, p, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: TDB
        real(wp),dimension(3),intent(in) :: p
        real(wp),dimension(3),intent(out) :: f

        real(wp),parameter :: pi = 3.141592653589793_wp
        real(wp),parameter :: E = 1358.0_wp
        real(wp),parameter :: c = 299792458.0_wp
        real(wp),parameter :: m = 180.0_wp
        real(wp),parameter :: d = 30.48_wp
        real(wp),parameter :: r = 0.88_wp
        real(wp),parameter :: s = 0.94_wp
        real(wp),parameter :: bf = 0.79_wp
        real(wp),parameter :: a0 = -(E/c)*(pi*(d/2)**2)/m
        real(wp),parameter :: ar = ((1+r*s)/2+(2*bf*r*(1-s))/3)*a0
        real(wp),parameter :: an = ar/(3986004.415D8/6378136.3_wp**2)
        real(wp),parameter :: re2 = (6378136.3_wp/149597870700.0_wp)**2
        real(wp),parameter :: AU = 149597870700.0_wp/6378136.3_wp
        real(wp),dimension(6) :: pves
        real(wp),dimension(3) :: pes
        real(wp),dimension(3) :: ps
        real(wp) :: des2
        real(wp) :: ds2

        call PLEPH(TDB, 11, 3, pves)

        pes = pves(1:3)
        ps = pes-p/AU
        des2 = sum(pes**2)
        ds2 = sum(ps**2)
        if ((sum((pes/sqrt(des2))*(ps/sqrt(ds2)))**2+re2/des2) > 1.0_wp) then
            if (sum(p*pes) > 0.0_wp) then
                f = an*pes/sqrt(des2)
            else
                f = 0.0_wp
            end if
        else
            f = an*pes/sqrt(des2)
        end if

    end subroutine srp1
