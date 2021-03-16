
    subroutine getfp(TDB, p, pvem, pves, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: TDB
        real(wp),dimension(3),intent(in) :: p
        real(wp),dimension(6),intent(in) :: pvem
        real(wp),dimension(6),intent(in) :: pves
        real(wp),dimension(3),intent(out) :: f

        real(wp),parameter :: AU = 149597870700.0_wp/6378136.3_wp
        real(wp),parameter :: GMm = 6.67430D-11*7.3477D22/3986004.415D8
        real(wp),parameter :: GMs = 6.67430D-11*1.9891D30/3986004.415D8
        real(wp),dimension(3) :: pem
        real(wp),dimension(3) :: pes
        real(wp),dimension(3) :: pm
        real(wp),dimension(3) :: ps
        real(wp),dimension(3) :: fm
        real(wp),dimension(3) :: fs

        pem = pvem(1:3)*AU
        pes = pves(1:3)*AU

        pm = pem-p
        ps = pes-p

        fm = GMm*(pm/sum(pm**2.0_wp)**1.5_wp-pem/sum(pem**2.0_wp)**1.5_wp)
        fs = GMs*(ps/sum(ps**2.0_wp)**1.5_wp-pes/sum(pes**2.0_wp)**1.5_wp)

        f = fm+fs

    end subroutine getfp
