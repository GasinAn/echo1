
    subroutine fcn(me,t,pv,vf)

    use dop853_module
    use dop853_constants
    
    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: pv
    real(wp),dimension(:),intent(out) :: vf

    integer,parameter :: n_max = 100
    real(wp),dimension(0:n_max,0:n_max) :: c_matrix, s_matrix
    real(wp),dimension(1:n_max,0:n_max-1) :: a
    real(wp),dimension(1:n_max) :: b
    real(wp),dimension(2:n_max,0:n_max-2) :: c
    real(wp),dimension(2:n_max) :: d
    real(wp),dimension(2:n_max,0:n_max-1) :: e
    real(wp),dimension(0:n_max,0:n_max) :: p
    real(wp),dimension(3) :: fe
    real(wp),dimension(3) :: fp

    call readegm(n_max, c_matrix, s_matrix)
    call getabcdep(n_max, a, b, c, d, e, p)
    call getfe(t, pv(1), pv(2), pv(3), &
               n_max, a, b, c, d, e, c_matrix, s_matrix, p, fe)
    call getfp(2400000.5D0+t, pv(1:3), fp)

    vf(1:3) = pv(4:6)
    vf(4:6) = fe+fp

    end subroutine fcn
