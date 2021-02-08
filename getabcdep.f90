
    subroutine getabcdep(n_max, a, b, c, d, e, p)

        use iso_fortran_env,    only: wp => real64

        implicit none

        integer,intent(in) :: n_max
        real(wp),dimension(1:n_max,0:n_max-1),intent(out) :: a
        real(wp),dimension(1:n_max),intent(out) :: b
        real(wp),dimension(2:n_max,0:n_max-2),intent(out) :: c
        real(wp),dimension(2:n_max),intent(out) :: d
        real(wp),dimension(2:n_max,0:n_max-1),intent(out) :: e
        real(wp),dimension(0:n_max,0:n_max),intent(out) :: p

        integer :: n, m
        real(wp) :: n_wp, m_wp

        do n = 1, n_max
            do m = 0, n-1
                n_wp = n
                m_wp = m
                a(n,m) = sqrt(((2*n_wp+1)*(2*n_wp-1))/((n_wp+m_wp)*(n_wp-m_wp)))
            end do
        end do

        b(1) = sqrt(3.0_wp)
        do n = 2, n_max
            n_wp = n
            b(n) = sqrt((2*n_wp+1)/(2*n_wp))
        end do

        do n = 2, n_max
            do m = 0, n-2
                c(n,m) = -a(n,m)/a(n-1,m)
            end do
        end do

        do n = 2, n_max
            d(n) = -(n+1)
        end do

        do n = 2, n_max
            do m = 0, n-1
                e(n,m) = -(2*n+1)/a(n,m)
            end do
        end do

        p = 0.0_wp
        p(0,0) = 1.0_wp

    end subroutine getabcdep
