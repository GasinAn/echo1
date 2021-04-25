
    subroutine pv2m(p, v)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),dimension(3),intent(in) :: p, v

        real(wp),parameter :: tau = 2*3.141592653589793_wp
        real(wp) :: a, e
        real(wp) :: px, py
        real(wp) :: cosEa, sinEa, Ea, M
        real(wp),dimension(3) :: vh, ve
        real(wp),dimension(3) :: vhbar, vebar, vpbar

        vh = [p(2)*v(3)-p(3)*v(2), &
              p(3)*v(1)-p(1)*v(3), &
              p(1)*v(2)-p(2)*v(1)]
        ve = [(v(2)*vh(3)-v(3)*vh(2))-p(1)/sqrt(sum(p**2)), &
              (v(3)*vh(1)-v(1)*vh(3))-p(2)/sqrt(sum(p**2)), &
              (v(1)*vh(2)-v(2)*vh(1))-p(3)/sqrt(sum(p**2))]

        a = sum(vh**2)/(1-sum(ve**2))
        e = sqrt(sum(ve**2))

        vhbar = vh/sqrt(sum(vh**2))
        vebar = ve/sqrt(sum(ve**2))
        vpbar = [vhbar(2)*vebar(3)-vhbar(3)*vebar(2), &
                 vhbar(3)*vebar(1)-vhbar(1)*vebar(3), &
                 vhbar(1)*vebar(2)-vhbar(2)*vebar(1)]

        px = sum(p*vebar)
        py = sum(p*vpbar)
        cosEa = px/a+e
        sinEa = py/a/sqrt(1-sum(ve*ve))
        Ea = atan2(sinEa,cosEa)
        M = Ea-e*sinEa
        print *, M/tau+1

    end subroutine pv2m
