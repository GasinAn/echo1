﻿
    subroutine getfe(t, p_gcrs, n_max, a, b, c, d, e, p, c_m, s_m, f)

        use iso_fortran_env,    only: wp => real64

        implicit none

        integer,intent(in) :: n_max
        real(wp),intent(in) :: t
        real(wp),dimension(3),intent(in) :: p_gcrs
        real(wp),dimension(1:n_max,0:n_max-1),intent(in) :: a
        real(wp),dimension(1:n_max),intent(in) :: b
        real(wp),dimension(2:n_max,0:n_max-2),intent(in) :: c
        real(wp),dimension(2:n_max),intent(in) :: d
        real(wp),dimension(2:n_max,0:n_max-1),intent(in) :: e
        real(wp),dimension(0:n_max,0:n_max),intent(inout) :: p
        real(wp),dimension(0:n_max,0:n_max),intent(in) :: c_m, s_m
        real(wp),dimension(3),intent(out) :: f

        integer :: n
        real(wp),dimension(3,3) :: rc2t, rt2c
        real(wp),dimension(3) :: p_itrs
        real(wp) :: r, cth, sth, cph, sph
        real(wp),dimension(0:n_max) :: mph, cmph, smph
        real(wp),dimension(0:n_max) :: ccss, sccs, fr, fth, fph
        real(wp),dimension(3) :: df

        call getrc2trt2c(t, rc2t, rt2c)

        p_itrs = [sum(rc2t(1,:)*p_gcrs), &
                  sum(rc2t(2,:)*p_gcrs), &
                  sum(rc2t(3,:)*p_gcrs)]

        r = sqrt(sum(p_itrs**2))
        cth = p_itrs(3)/r
        sth = sqrt(sum(p_itrs(1:2)**2))/r

        p(1,0:1) = [a(1,0)*cth*p(0,0),b(1)*sth*p(0,0)]
        do n = 2, n_max
            p(n,0:n) = [a(n,0:n-2)*cth*p(n-1,0:n-2)+c(n,0:n-2)*p(n-2,0:n-2), &
                        a(n,n-1)*cth*p(n-1,n-1), &
                        b(n)*sth*p(n-1,n-1)]
        end do

        mph = [0:n_max]*atan2(p_itrs(2),p_itrs(1))
        cmph = cos(mph)
        smph = sin(mph)
        cph = cmph(1)
        sph = smph(1)

        f = -p_itrs/r**3
        do n = 2, n_max
            ccss(0:n) = c_m(n,0:n)*cmph(0:n)+s_m(n,0:n)*smph(0:n)
            sccs(0:n) = s_m(n,0:n)*cmph(0:n)-c_m(n,0:n)*smph(0:n)
            fr = d(n)*p(n,0:n)*ccss(0:n)
            fth = [(n*cth*p(n,0:n-1)+e(n,0:n-1)*p(n-1,0:n-1))/sth*ccss(0:n-1), &
                   n*cth*p(n,n)/sth*ccss(n)]
            fph = [0:n]*p(n,0:n)/sth*sccs(0:n)
            df = [sum(fr(0:n)), sum(fth(0:n)), sum(fph(0:n))]
            df = [(df(1)*sth+df(2)*cth)*cph-df(3)*sph, &
                  (df(1)*sth+df(2)*cth)*sph+df(3)*cph, &
                  (df(1)*cth-df(2)*sth)]
            f = f+df/r**(n+2)
        end do

        f = [sum(rt2c(1,:)*f), &
             sum(rt2c(2,:)*f), &
             sum(rt2c(3,:)*f)]

    end subroutine getfe
