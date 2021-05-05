!*****************************************************************************************
!> author: GasinAn
!
!  Simplified Modern Fortran Edition of the DOP853 ODE Solver.
!
!### License
!
!        Simplified Modern Fortran Edition of the DOP853 ODE Solver
!        https://github.com/GasinAn/easydop853
!
!        Copyright (c) 2020, GasinAn
!        All rights reserved.
!
!        Redistribution and use in source and binary forms, with or without modification,
!        are permitted provided that the following conditions are met:
!
!        * Redistributions of source code must retain the above copyright notice, this
!          list of conditions and the following disclaimer.
!
!        * Redistributions in binary form must reproduce the above copyright notice, this
!          list of conditions and the following disclaimer in the documentation and/or
!          other materials provided with the distribution.
!
!        * The names of its contributors may not be used to endorse or promote products
!          derived from this software without specific prior written permission.
!
!        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
!        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
!        ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!*****************************************************************************************

    module easydop853_module

    use dop853_module
    use dop853_constants

    implicit none

    contains
    
    subroutine easydop853(fcn,x,xf,y)

    implicit none

    procedure(deriv_func)               :: fcn
    !! subroutine computing the value of \(dy/dx=f(x,y)\)
    real(wp),intent(inout)              :: x
    !! `x` value (input is initial value and output is final value)
    real(wp),intent(in)                 :: xf
    !! endpoint of integration (final value of `x`)
    real(wp),dimension(:),intent(inout) :: y
    !! `y` value (input is initial value and output is final value)

    integer               :: nstiff = 1
    !! nstiff parameter for stiffness detection,
    !! which will occur at step 1*nstiff, 2*nstiff, 3*nstiff ... if nstiff>0
    !! and will not occur if nstiff<=0
    integer               :: nmax   = 2250000
    !! maximal number of allowed steps
    !real(wp),dimension(1) :: rtol   = 1.0e-10_wp
    !real(wp),dimension(1) :: rtol   = 2.0e-12_wp
    real(wp),dimension(1) :: rtol   = 1.0e-12_wp
    !! relative tolerance
    !real(wp),dimension(1) :: atol   = 1.0e-24_wp
    real(wp),dimension(1) :: atol   = 0.0_wp
    !! absolute tolerance

    type(dop853_class) :: prop
    logical :: status_ok
    integer :: idid

    call prop%initialize(n=size(y),fcn=fcn,nstiff=nstiff,nmax=nmax,&
                         status_ok=status_ok)
    if (.not.status_ok) error stop 'initialization error'

    call prop%integrate(x,y,xf,rtol,atol,iout=0,idid=idid)
    if (idid<0) error stop 'integration failure'

    end subroutine easydop853

    end module easydop853_module
