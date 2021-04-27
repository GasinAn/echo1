
    subroutine srp2(TDB, p, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: TDB
        real(wp),dimension(3),intent(in) :: p
        real(wp),dimension(3),intent(out) :: f

        real(wp),parameter :: pi = 3.141592653589793_wp
        real(wp),parameter :: E = 1358.0_wp
        !real(wp),parameter :: E = 1350.0_wp
        real(wp),parameter :: c = 299792458.0_wp
        !real(wp),parameter :: m = 180.0_wp
        real(wp),parameter :: m = 156.995_wp*0.45359237_wp
        !real(wp),parameter :: m = (156.995_wp-33.34_wp)*0.45359237_wp
        real(wp),parameter :: d = 30.48_wp
        real(wp),parameter :: r = 0.88_wp
        real(wp),parameter :: s = 0.94_wp
        real(wp),parameter :: bf = 0.79_wp
        real(wp),parameter :: bb = 0.55_wp
        real(wp),parameter :: ef = 0.05_wp
        real(wp),parameter :: eb = 0.05_wp
        real(wp),parameter :: a0 = -(E/c)*(pi*(d/2)**2)/m
        real(wp),parameter :: k = (1+r*s)/2+bf*(1-s)*r*2/3
        real(wp),parameter :: kp = k+(1-r)*(ef*bf-eb*bb)/(ef+eb)*2/3
        real(wp),parameter :: an = k*a0/(3986004.415D8/6378136.3_wp**2)
        !real(wp),parameter :: an = kp*a0/(3986004.415D8/6378136.3_wp**2)
        real(wp),parameter :: AU = 149597870700.0_wp/6378136.3_wp

        real(wp),dimension(6) :: pves
        real(wp),dimension(6) :: pvem
        real(wp),dimension(3) :: pes
        real(wp),dimension(3) :: pms
        real(wp),dimension(3) :: pm
        real(wp) :: dpes
        real(wp) :: dp
        real(wp) :: dpms
        real(wp) :: dpm
        real(wp) :: th1
        real(wp) :: th2
        real(wp) :: th3
        real(wp) :: th4
        real(wp) :: costh12
        real(wp) :: costh34
        real(wp),dimension(3) :: ps

        call PLEPH(TDB, 11, 3, pves)
        call PLEPH(TDB, 10, 3, pvem)

        !pes = pves(1:3)
        !ps = pes-p/AU
        !des2 = sum(pes**2)
        !ds2 = sum(ps**2)
        !if ((sum((pes/sqrt(des2))*(ps/sqrt(ds2)))**2+re2/des2) > 1.0_wp) then
        !    if (sum(p*pes) > 0.0_wp) then
        !        f = an*pes/sqrt(des2)
        !    else
        !        f = 0.0_wp
        !    end if
        !else
        !    f = an*pes/sqrt(des2)
        !end if

        pes = pves(1:3)
        pms = pes-pvem(1:3)
        pm = p-pvem(1:3)*AU

        dpes = sqrt(sum(pes**2))
        dp = sqrt(sum(p**2))
        dpms = sqrt(sum(pms**2))
        dpm = sqrt(sum(pm**2))

        th1 = asin(6371000.0_wp/(dpes*149597870700.0_wp))
        !th1 = asin(6421000.0_wp/(dpes*149597870700.0_wp))
        th2 = acos((6371000.0_wp/6378136.3_wp)/dp)
        !th2 = acos((6421000.0_wp/6378136.3_wp)/dp)
        th3 = asin(1737400.0_wp/(dpms*149597870700.0_wp))
        th4 = acos((1737400.0_wp/6378136.3_wp)/dpm)

        costh12 = sum((p/dp)*(pes/dpes))
        costh34 = sum((pm/dpm)*(pms/dpms))

        if ((costh12>sin(th1-th2)).and.(costh34>sin(th3-th4))) then
            !f = an*(pes/dpes)
            ps = pes-p/AU
            f = an*(ps/norm2(ps))
        else
            f = 0
        end if

    end subroutine srp2
