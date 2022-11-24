common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8
common/cosmo_today/om_m0, om_v0
common/norm/pfac
common/pdat/x_pdat(1000), y_pdat(1000), n_pdat
common/sigdat/datlook_nu(1000), datlook_m(1000)
real*8 sum_check, sum_bar, tot, tot_ps

open (6, file = 'C:\\Users\\lexgu\\Desktop\\MPhys\\Data\\output.txt', status = 'replace')

pi = 3.141592654

! Set parameters
om_m0 = 0.3
om_v0 = 0.7
sig80 = 0.8
p_index = 0.96

h = 0.7
om_b = 0.045

rho0 = 2.78e11 * om_m0

! Comoving mass density. No powers of h: absorbed in length and mass units
do iz = 0, 11
      z = 1.0 * iz

      om_m = om_m0
      om_v = om_v0
      sig8 = sig80/grow(z)

      write(6,*) 'z = ', z, ' sigma_8(z) = ', sig8

      om_m1 = omega_m(z)
      om_v1 = omega_v(z)
      om_m = om_m1
      om_v = om_v1

      pfac = 1.

      ! need pfac=1 before iterating normalization to match sigma_8
      pfac = sig8/sigint(8.)
      pfac = pfac * pfac

      ! Initialize nu-mass look-up table
      call init_look()

      ! Loop backwards through look-up table to get integral density
      tot = 0.d0
      tot_ps = 0.d0

      sum_check = 0.d0
      sum_bar = 0.d0

      do i = 999, 1, -1
            rnu1 = datlook_nu(i)
            rnu2 = datlook_nu(i + 1)
            rm1 = datlook_m(i)
            rm2 = datlook_m(i + 1)
            rnu = (rnu1 + rnu2)/2.

            rm = (rm1 + rm2)/2./h

            ! Here we convert mass from internal units to M_sun
            df = fcoll(rnu1) - fcoll(rnu2)

            sum_check = sum_check+df
            sum_bar = sum_bar + df * brat(log10(rm))

            fmult = df/log(rm2/rm1)
            fac = abs(log(rnu1/rnu2)/log(rm2/rm1))
            fmult_ps = fac * sqrt(2./pi) * rnu * exp(-rnu * rnu/2.)

            den = fmult * (rho0/rm) * log(10.)
            den_ps = fmult_ps * (rho0/rm) * log(10.)

            tot = tot + fmult * (rho0/rm) * log(rm2/rm1)
            tot_ps = tot_ps + fmult_ps * (rho0/rm) * log(rm2/rm1)
      
            write(6,*) i, rm, rnu, fmult, fcoll(rnu), tot
            ! write(6,*) sum_check, sum_bar
      enddo
enddo

close(6)

end

function p_cdm(rk)
      common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8

      p_cdm = p_full(rk)

      return
end

subroutine init_look()
      common/cosmo_today/om_m0, om_v0
      common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8
      common/sigdat/datlook_nu(1000), datlook_m(1000)

      pi = 3.141592654

      ! cover a log grid in filter radius - units comoving Mpc/h
      do i = 1, 1000
            rfil = -2. + 4. * (i - 1)/999.
            rfil = 10. **rfil
            rmass = rfil **3 * 4 * pi * 2.78e11 * om_m0/3.
            rnu = 1.686/sigint(rfil)

            datlook_nu(i) = rnu
            datlook_m(i) = rmass
      enddo

      return
end

function sigint(r)
      implicit real*8 (a - h, o - z)
      real*4 sigint, r, p_cdm, rk4

      nint = 899
      sum1 = 0.d0

      do i = 1, nint
            t = (float(i) - 0.5)/float(nint)
            y = -1.d0 + 1.d0/t

            d2 = p_cdm(rk4)

            x = y * r
            w = (3./x/x/x) * (sin(x) - x * cos(x))

            sum1 = sum1 + w * w * d2/y/t/t
      enddo

      sum1 = sum1/float(nint)
      sigint = sqrt(sum1)

      return
end

function grow(z)
      common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8

      onow_m = omega_m(z)
      onow_v = omega_v(z)

      grow = gg(onow_m, onow_v)
      grow0 = gg(om_m, om_v)

      ! grow=(1 + z)*grow0/grow
      grow = grow0/grow

      return
end

function omega_m(z)
      common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8

      aa = 1/(1 + z)
      omega_t = 1.0 + (om_m + om_v - 1.0)/(1 - om_m - om_v + om_v * aa * aa + om_m/aa)
      omega_m = omega_t * om_m/(om_m + om_v * aa * aa * aa)

      return
end
            
function omega_v(z)
      common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8

      aa = 1/(1 + z)
      omega_t = 1.0 + (om_m + om_v - 1.0)/(1 - om_m - om_v + om_v * aa * aa + om_m/aa)
      omega_v = omega_t * om_v/(om_v + om_m/aa/aa/aa)

      return
end

function gg(om_m, om_v)
      ! gg = (5./2.) * om_m/(om_m **(4./7.) - om_v + (1 + om_m/2.) * (1 + om_v/70.))
      x = om_v **(1./3.)
      gg = x * (1 - x **(1.91)) **(0.82) + 1.437 * (1 - (1 - x **3) **(2./3.))

      ! For flat case only, and now gg is entire growth, not growth/a
      return
end

function p_full(rk)
      common/cospar/om_m, om_v, om_b, h, p_index, gams, sig8
      common/norm/pfac

      ! approximate COBE normalization allowing for p_index
      tilt = p_index-1

      if(om_v.gt.0.) then
            dh = 1.94e-5 * (om_m **(-0.785 - 0.05 * log(om_m))) * exp(-0.95 * tilt - 0.169 * tilt * tilt)
      else
            dh = 1.95e-5 * (om_m **(-0.35 - 0.19 * log(om_m) - 0.17 * tilt)) * exp(-tilt - 0.14 * tilt * tilt)
      endif

      p_full = dh * dh * (3000. * rk) **(3 + p_index) * tk_eh(rk) * tk_eh(rk) * pfac

      return
end

function tk_eh(yy)
      common/cospar/om_mt, om_vt, om_b, h, p_index, gams, sig8
      common/cosmo_today/om_m, om_v

      ! the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
      ! remember I use k/h, whereas they use pure k
      
      ! om_m is the total matter density parameter - i.e. CDM + baryons
      rk = yy * h

      e = exp(1.)
      
      thet = 2.728/2.7
      b1 = 0.313 * (om_m * h * h) **(-0.419) * (1 + 0.607 * (om_m * h * h) **0.674)
      b2 = 0.238 * (om_m * h * h) **0.223
      zd = 1291 * (1 + b1 * (om_b * h * h) **b2) * (om_m * h * h) **0.251/(1 + 0.659 * (om_m * h * h) **0.828)
      ze = 2.50e4 * om_m * h * h/thet **4
      rd = 31500 * om_b * h * h/thet **4/zd
      re = 31500 * om_b * h * h/thet **4/ze
      rke = 7.46e-2 * om_m * h * h/thet **2
      s = (2./3./rke) * sqrt(6./re) * log((sqrt(1 + rd) + sqrt(rd + re))/(1 + sqrt(re)))
      rks = 1.6 * ((om_b * h * h) **0.52) * ((om_m * h * h) **0.73) * (1 + (10.4 * om_m * h * h) **(-0.95))

      q = rk/13.41/rke
      
      y = (1 + ze)/(1 + zd)
      g = y * (-6 * sqrt(1 + y)+(2 + 3 * y) * log((sqrt(1 + y) + 1)/(sqrt(1 + y) - 1)))
      ab = g * 2.07 * rke * s/(1 + rd) **(0.75)

      a1 = (46.9 * om_m * h * h) **0.670 * (1 + (32.1 * om_m * h * h) **(-0.532))
      a2 = (12.0 * om_m * h * h) **0.424 * (1 + (45.0 * om_m * h * h) **(-0.582))
      ac = (a1 **(-om_b/om_m)) * (a2 **(-(om_b/om_m) **3))

      b1 = 0.944/(1 + (458 * om_m * h * h) **(-0.708))
      b2 = (0.395 * om_m * h * h) **(-0.0266)
      bc = 1/(1 + b1 *((1 - om_b/om_m) **b2-1))

      f = 1/(1 + (rk * s/5.4) **4)

      c1 = 14.2 + 386/(1 + 69.9 * q **1.08)
      c2 = 14.2/ac + 386/(1 + 69.9 * q **1.08)
      tc = f * log(e + 1.8 * bc * q)/(log(e + 1.8 * bc * q) + c1 * q * q) &
            + (1 - f) * log(e + 1.8 * bc * q)/(log(e + 1.8 * bc * q) + c2 * q * q)

      bb = 0.5 + (om_b/om_m) + (3 - 2 * om_b/om_m) * sqrt((17.2 * om_m * h * h) **2 + 1)
      bn = 8.41 * (om_m *h *h) **0.435
      ss = s/(1 + (bn/rk/s) **3) **(1./3.)
      tb = log(e + 1.8 * q)/(log(e + 1.8 * q) + c1 * q * q)/(1 + (rk * s/5.2) **2)
      tb = (tb + ab * exp(-(rk/rks) **1.4)/(1 + (bb/rk/s) **3)) * sin(rk * ss)/rk/ss

      tk_eh = (om_b/om_m) * tb + (1 - om_b/om_m) * tc

      return
end

function fcoll(x)
      ! JAP 2007 expression for integral collapse fraction
      fcoll = exp(-0.412 * x * x)/(1 + 1.529 * x **0.704)

      return
end

function brat(x)
      y = x - 12.9
      brat = 1./(1. + exp(-y/1.12))

      return
end
