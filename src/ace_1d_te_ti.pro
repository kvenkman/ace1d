PRO ace_1d_te_ti, model, pconst, zion, zmajnow, zionnow, eheating, cooltermse, heattermsi, cooltermsi, heatterms

; Refer to TIE-GCM documentation Eq. 6.54 onwards
; Boundary conditions - 
; For LBC, Te and Ti are tied to the neutral temperature
; For UBC, we use d/dz(dTi/dz) = 0, and specify a downward heatflux for Te (Smithtro dissertation, pg. 12-13)
;   heat_flux = 4.4e9 * 1.6e-12; ergs cm^-2 s^-1
    heat_flux = 3e9 * 1.6e-12; ergs cm^-2 s^-1  ; Roble [1987]
    sc_h = pconst.boltz*zmajnow.tn/(pconst.grav*zmajnow.barm/pconst.avo) ; Mean scale height in cm
	
;   i = 90.*(!pi/180.) ; Magnetic dip angle., Roble [1987] 
    i = 75.*(!pi/180.) ; Magnetic dip angle., Smithtro [2005]
    sin_i_sq = (sin(i))^2.
    
    ; Calculating electron thermal conductivity [Rees and Roble, 1975; Eq. 40]
    ; Also, read paragraph beneath Eq. 5.146 of Ionospheres, Schunk & Nagy (2009)
    ; for a funny insight into why this equation is "wrong".
    q1 = 2.82e-17*(zionnow.te^0.5) - 3.41e-21*(zionnow.te^1.5)
    q2 = 2.2e-16 + 7.92e-18*(zionnow.te^0.5)
    q3 = 3.4e-16
    qs = q1*zmajnow.n2den + q2*zmajnow.o2den + q3*zmajnow.oden    
    lambda = 1.6e-12 * 7.5e5/(1. + 3.22e4*qs*(zionnow.te^2./zionnow.e)) ; ergs cm^-1 s^-1 K ^-1

    t1 = sin_i_sq/sc_h
    t2 = (2./7.)*(lambda/sc_h)
    t3 = -(cooltermse.neutrals_implicit + cooltermse.ions_implicit)/zionnow.te^2.5

    dt2 = deriv(model.zp, t2)
    t4 = t1*dt2
    t5 = t1*t2

;;;; Electron Temperature solver
      
    p = fltarr(model.nlev)
    q = fltarr(model.nlev)
    r = fltarr(model.nlev)
    rhs = fltarr(model.nlev)
    
    p = t5/(model.dz^2.) - t4/(2.*model.dz)
    q = t3 - 2.*t5/(model.dz^2.)
    r = t5/(model.dz^2.) + t4/(2.*model.dz)
    rhs = (-1.)*(eheating + cooltermse.neutrals_implicit*zmajnow.tn + cooltermse.ions_implicit*zionnow.ti) + cooltermse.explicit

    ; Boundary conditions
    ; LBC
    p[0] = 0.
    q[0] = 1.
    r[0] = 0.
    rhs[0] = zmajnow.tn[0]^3.5
    
    ; UBC
    p[model.nlev-1] = p[model.nlev-1] + r[model.nlev-1]
    rhs[model.nlev-1] = rhs[model.nlev-1] - (2.*model.dz)*heat_flux*r[model.nlev-1]*(7./2.)*(sc_h[model.nlev-1]/lambda[model.nlev-1]) ; zionnow.te[model.nlev-1]; 
    r[model.nlev-1] = 0.    
    
    te_update = (trisol(p, q, r, rhs))^(2./7.) > zmajnow.tn
;    te_update = te_update^(2./7.)
    zion.te = te_update

;;;; Ion Temperature solver 
;    p = fltarr(model.nlev)
;    q = fltarr(model.nlev)
;    r = fltarr(model.nlev)
;    rhs = fltarr(model.nlev)
    
;    pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn)
;    rho = pkt*(zmajnow.barm/pconst.avo)

; Ion temperature is set assuming a simple thermodynamic balance
    ti_update = (zionnow.te*heattermsi.ei + zmajnow.tn*cooltermsi.neut)/(cooltermsi.neut + heattermsi.ei) > 1e-24
    ti_update = ti_update > zmajnow.tn
    zion.ti = ti_update
    
; cgplot, zionnow.te, zmajnow.zz
; cgoplot, zion.te, zmajnow.zz, linestyle = 2
; cgoplot, zionnow.ti, zmajnow.zz, color = 'blue'
; cgoplot, zion.ti, zmajnow.zz, color = 'blue', linestyle = 2
 
; cgoplot, zmajnow.tn, zmajnow.zz, color = 'red'

; wait, 0.1
END
