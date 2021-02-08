
    subroutine readegm(n_max,c_matrix,s_matrix)

        use iso_fortran_env,    only: wp => real64

        implicit none

        integer,intent(in) :: n_max
        real(wp),dimension(0:n_max,0:n_max),intent(out) :: c_matrix, s_matrix

        integer :: i, j

        integer :: n, m
        real(wp) :: c, s, dc, ds

        c_matrix = 0.0_wp
        s_matrix = 0.0_wp

        open(10,file='EGM2008_to2190_ZeroTide')

        do i = 1, (n_max+4)*(n_max-1)/2
            read(10,'(2i5,2d25.15,2d20.10)') n, m, c, s, dc, ds
            c_matrix(n,m) = c
            s_matrix(n,m) = s
        end do

        close(10)

    end subroutine readegm
