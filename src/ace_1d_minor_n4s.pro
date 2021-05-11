FUNCTION ace_1d_minor_n4s, zminornow, zmajnow, prod, loss, model, mass, pconst
; Solving for minor species mass mixing ratios by accounting for 
; production, loss and diffusion

; The lower boundary condition is defined as PCE. 
; Certain models define a flux for N(4S) as the LBC, but we do not presently do so.
; The upper boundary condition is defined to be diffusive equilibrium
    
    k_e = 5e-6*exp(-7-model.zp)*exp(-model.zp)
    p0  = model.p0 ; 5e-4, in dyne/cm^2
    p00 = 1e6
    t00 = 273.
	h0 = pconst.boltz*t00/(mass.n2*pconst.grav/pconst.avo) ; characteristic scale height of N2        
	d0  = 0.2
    tau = (1. + fltarr(model.nlev))*p0*(h0^2.)/(p00*d0)
    t0 = 273.
    
    phi12 = 1.35
    phi13 = 1.11
    phi21 = 0.673
    phi23 = 0.769
    
    phix1 = 0.651
    phix2 = 0.731
    phix3 = 0.741
    
    pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn)
    prod_mmr = prod*mass.n4s/(zmajnow.barm*pkt) ; Converting production rate from molecules cm^-3 s^-1 to mass mixing ratio
    loss_mmr = loss ; Loss frequency remains the same in mmr

    alphax1 = phix1 - phix3
    alphax2 = phix2 - phix3
    
    salfa12 = phi12 - phi13
    salfa21 = phi21 - phi23
    salfax1 = phix1 - phix3
    salfax2 = phix2 - phix3

    alpha11 = - (phi13 + salfa12*zmajnow.o)
    alpha12 = salfa12*zmajnow.o2
    alpha21 = salfa21*zmajnow.o
    alpha22 = - (phi23 + salfa21*zmajnow.o2)
    
    dbarm = deriv(model.zp, zmajnow.barm)
    do2   = deriv(model.zp, zmajnow.o2)
    do1   = deriv(model.zp, zmajnow.o )
    
    a_m = zmajnow.barm/(mass.n2)*(t0/zmajnow.tn)^0.25
    a_m = a_m/(phix3 + salfax1*zmajnow.o2 + salfax2*zmajnow.o)
    
    e_x =  ((salfax1*alpha22 - salfax2*alpha21)*(do2 - (1. - (mass.o2 + dbarm)/zmajnow.barm)*zmajnow.o2) + $
           (salfax2*alpha11 - salfax1*alpha12)*(do1 - (1. - (mass.o  + dbarm)/zmajnow.barm)*zmajnow.o ))  $
           /(alpha11*alpha22 - alpha12*alpha21) + $
           (1. - mass.n4s/zmajnow.barm - dbarm/zmajnow.barm)
           
    da_m = deriv(model.zp, a_m)
    damex = deriv(model.zp, a_m*e_x)
    dk_e = deriv(model.zp, k_e)
    dkebarm_m = deriv(model.zp, (k_e/zmajnow.barm)*dbarm)
    
    t1 = ((1./tau)*a_m + k_e)
    t2 = ((1./tau)*da_m - (1./tau)*a_m*e_x + dk_e + k_e*dbarm/zmajnow.barm)
    t3 = (-(1./tau)*damex - loss*exp(-model.zp) + dkebarm_m)
    
    p   = fltarr(model.nlev)
    q   = fltarr(model.nlev)
    r   = fltarr(model.nlev)
    rhs = fltarr(model.nlev)
    
    p = (t1/(model.dz^2.) - t2/(2.*model.dz))
    q = (-2.*t1/(model.dz^2.) + t3 - exp(-model.zp)/model.del_time)
    r = (t1/(model.dz^2.) + t2/(2.*model.dz))
    rhs = -(exp(-model.zp)*zminornow.n4s_mmr/model.del_time + prod_mmr*exp(-model.zp))

    ; LBC, PCE
    p[0] = 0.
    q[0] = 1.
    r[0] = 0.
    rhs[0] = prod_mmr[0]/loss_mmr[0]
;	flux = -1e-4
	
;	q[0] = q[0] - p[0]*(2.*model.dz)*e_x[0]
;	r[0] = r[0] + p[0]
;	rhs[0] = rhs[0] + p[0]*flux*(2.*model.dz)/a_m[0]
;	p[0] = 0.
	
    ; UBC
    p[model.nlev-1] = p[model.nlev-1] + r[model.nlev-1]
    q[model.nlev-1] = q[model.nlev-1] + (2.*model.dz)*e_x[model.nlev-1]*r[model.nlev-1]
    r[model.nlev-1] = 0.
    
    n4s_mmr_upd = trisol(p, q, r, rhs) > 0.
    
    return, n4s_mmr_upd
    END
