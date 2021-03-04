
    subroutine pv2ae(p, v, a, e)

        use iso_fortran_env,    only: wp => real64

        implicit none

        real(wp),dimension(3),intent(in) :: p, v
        real(wp),intent(out) :: a, e

        real(wp),dimension(3) :: vh, ve

        vh = [p(2)*v(3)-p(3)*v(2), &
              p(3)*v(1)-p(1)*v(3), &
              p(1)*v(2)-p(2)*v(1)]
        ve = [(v(2)*vh(3)-v(3)*vh(2))-p(1)/sqrt(sum(p**2)), &
              (v(3)*vh(1)-v(1)*vh(3))-p(2)/sqrt(sum(p**2)), &
              (v(1)*vh(2)-v(2)*vh(1))-p(3)/sqrt(sum(p**2))]
        a = sum(vh**2)/(1-sum(ve**2))
        e = sqrt(sum(ve**2))

    end subroutine pv2ae
