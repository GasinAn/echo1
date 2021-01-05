
    subroutine fcn(me,x,y,f)

    !! subroutine computing the value of \(dy/dx=f(x,y)\)

    use dop853_module
    use dop853_constants
    
    implicit none

    class(dop853_class),intent(inout) :: me
    !! dop853_class object
    real(wp),intent(in)               :: x
    !! independent variable \(x\)
    real(wp),dimension(:),intent(in)  :: y
    !! state vector \(y(x)\)
    real(wp),dimension(:),intent(out) :: f
    !! derivative vector \(f(x,y)=dy/dx\)

    f(1:3)=y(4:6)
    f(4:6)=-y(1:3)/sum(y(1:3)**2.0_wp)**1.5_wp

    end subroutine fcn
