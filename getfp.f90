
    subroutine getfp(TDB, p, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),intent(in) :: TDB
        real(wp),dimension(3),intent(in) :: p
        real(wp),dimension(3),intent(out) :: f

        real(wp),parameter :: AU = 149597870700.0_wp
        real(wp),parameter :: GMm = 6.67430D-11*7.3477D22
        real(wp),parameter :: GMs = 6.67430D-11*1.9891D30
        real(wp),dimension(6) :: pvme
        real(wp),dimension(6) :: pvse
        real(wp),dimension(3) :: pme
        real(wp),dimension(3) :: pse
        real(wp),dimension(3) :: pm
        real(wp),dimension(3) :: ps
        real(wp),dimension(3) :: fm
        real(wp),dimension(3) :: fs

        call PLEPH(TDB, 3, 10, pvme)
        call PLEPH(TDB, 3, 11, pvse)

        pme = pvme(1:3)*AU
        pse = pvse(1:3)*AU

        pm = pme-p
        ps = pse-p

        fm = GMm*(pm/sum(pm**2.0_wp)**1.5_wp-pme/sum(pme**2.0_wp)**1.5_wp)
        fs = GMs*(ps/sum(ps**2.0_wp)**1.5_wp-pse/sum(pse**2.0_wp)**1.5_wp)

        f = fm+fs

    end subroutine getfp
