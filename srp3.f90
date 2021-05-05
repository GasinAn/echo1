
    subroutine srp3(TDB, p, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: TDB
        real(wp),dimension(3),intent(in) :: p
        real(wp),dimension(3),intent(out) :: f

        real(wp),parameter :: pi = 3.141592653589793_wp
        real(wp),parameter :: Es = 1358.0_wp
        !real(wp),parameter :: Es = 1350.0_wp
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
        real(wp),parameter :: alpha = 0.3_wp
        real(wp),parameter :: etae1 = (2*alpha)/(3*pi**2)
        real(wp),parameter :: etae2 = (1-alpha)/(4*pi)
        real(wp),parameter :: Re = 6371000.0_wp
        !real(wp),parameter :: Re = 6421000.0_wp
        real(wp),parameter :: Ae = pi*(Re/6378136.3_wp)**2
        real(wp),parameter :: a0 = (Es/c)*(pi*(d/2)**2)/m
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
        real(wp),dimension(3) :: np
        real(wp),dimension(3) :: npes
        real(wp),dimension(3) :: npm
        real(wp),dimension(3) :: npms
        real(wp) :: costh12
        real(wp) :: costh34
        real(wp) :: th12
        real(wp) :: sinth12
        real(wp) :: etae
        real(wp),dimension(3) :: ps

        call PLEPH(TDB, 11, 3, pves)
        call PLEPH(TDB, 10, 3, pvem)

        pes = pves(1:3)
        pms = pes-pvem(1:3)
        pm = p-pvem(1:3)*AU

        dpes = sqrt(sum(pes**2))
        dp = sqrt(sum(p**2))
        dpms = sqrt(sum(pms**2))
        dpm = sqrt(sum(pm**2))

        th1 = asin(Re/(dpes*149597870700.0_wp))
        th2 = acos((Re/6378136.3_wp)/dp)
        th3 = asin(1737400.0_wp/(dpms*149597870700.0_wp))
        th4 = acos((1737400.0_wp/6378136.3_wp)/dpm)

        np = p/dp
        npes = pes/dpes
        npm = pm/dpm
        npms = pms/dpms
        costh12 = sum(np*npes)
        costh34 = sum(npm*npms)

        if ((costh12>sin(th1-th2)).and.(costh34>sin(th3-th4))) then
            th12 = acos(costh12)
            sinth12 = sqrt(1-costh12**2)
            etae = (etae1*((pi-th12)*costh12+sinth12)+etae2)*(Ae/(sum(p**2)))
            ps = pes-p/AU
            f = an*(etae*np-ps/norm2(ps))
            !f = an*(etae*np-npes)
            !f = an*(-npes)
        else
            f = 0
        end if

    end subroutine srp3
