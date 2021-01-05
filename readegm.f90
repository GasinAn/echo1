
    subroutine readegm(n_max, nm, cs)

        use iso_fortran_env,    only: wp => real64

        implicit none

        integer,intent(in) :: n_max
        integer,dimension(2,(n_max+4)*(n_max-1)/2),intent(out) :: nm
        real(wp),dimension(2,(n_max+4)*(n_max-1)/2),intent(out) :: cs

        integer :: i, j

        integer :: n, m
        real(wp) :: c, s, dc, ds
        real(wp) :: a

        open(10,file='EGM2008_to2190_ZeroTide')

        do i = 1, (n_max+4)*(n_max-1)/2
            read(10,'(2i5,2d25.15,2d20.10)') n, m, c, s, dc, ds
            nm(:,i) = [n,m]
            cs(:,i) = [c,s]
            a = 2*n+1
            if (m.ne.0) then
                a = a*2.0_wp
            end if
            do j = n-m+1, n+m
                a = a/j
            end do
            cs(:,i) = cs(:,i)*sqrt(a)
        end do
        
        close(10)

    end subroutine readegm
