FUNCTION ace_1d_minor_solv, x, x_mass, zmajnow, prod, loss, k_m, flux_term, model, mass, pconst,i
; Solving for minor species mass mixing ratios by accounting for 
; production, loss and diffusion

; The lower boundary condition is defined as PCE. 
; Certain models define a flux for N(4S) as the LBC, but we do not presently do so.
; The upper boundary condition is defined to be diffusive equilibrium
    p   = fltarr(model.nlev)
    q   = fltarr(model.nlev)
    r   = fltarr(model.nlev)
    rhs = fltarr(model.nlev)
   
    h = pconst.boltz*zmajnow.tn/(pconst.grav*zmajnow.barm/pconst.avo) ; Scale height in cm
    pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn) ; Number density
    
    dtn = deriv(model.zp, zmajnow.tn)
    sc_h = pconst.boltz*zmajnow.tn/(pconst.grav*zmajnow.barm/pconst.avo)
    
    k_e = 5e-6*exp(-7-model.zp)*(sc_h^2.) ; Roble, 1987

    ; if(x_mass EQ 30) then k_e = 0.*k_e
    ;k_e = 1e6*exp(-7-model.zp) ; Eddy diffusion coefficient from Smithtro, 2005
    
    lbc = (x[0] + prod[0]*model.del_time)/(1.+loss[0]*model.del_time)
    
    t1 = (k_m + k_e)/h
    t2 = (1./h)*(k_m*(x_mass/zmajnow.barm + dtn/zmajnow.tn - 0.*h*flux_term) + k_e*(1. + dtn/zmajnow.tn))

    dt1 = deriv(model.zp, t1)
    dt2 = deriv(model.zp, t2)

    t3 = dt1 + t2
    
    p = t1/model.dz^2. - t3/(2.*model.dz)
    q = dt2 - 2.*t1/model.dz^2. - h/model.del_time - h*loss
    r = t1/model.dz^2. + t3/(2.*model.dz)

    rhs = (-1.*h)*(x/model.del_time + prod)

    ; Setting LBC, using PCE for N(4S) and specifying a downward flux for NO [Roble, 1987]
    if (x_mass EQ mass.n4s) then begin
        p[0] = 0.
        q[0] = 1.
        r[0] = 0.
        rhs[0] = lbc
    endif else begin
        ; LBC for NO is a small downward flux the allows for NO densities to peak at 105 km
        flux = 1e7
        q[0] = p[0]*(t2[0]/t1[0])*(2.*model.dz) + q[0]
        r[0] = p[0] + r[0]
        rhs[0] = rhs[0] + p[0]*flux*(2.*model.dz)/t1[0]

;        p[0] = 0.
;        q[0] = 1.
;        r[0] = 0.
;        rhs[0] = lbc; 1e8
    endelse
   
    ; UBC, diffusive equilibrium
    p[model.nlev-1] = p[model.nlev-1] + r[model.nlev-1]
    q[model.nlev-1] = q[model.nlev-1] - (2.*model.dz)*r[model.nlev-1]*(t2[model.nlev-1]/t1[model.nlev-1])
    r[model.nlev-1] = 0.
    x_np_cm3 = trisol(p, q, r, rhs) > 0.
    return, x_np_cm3 ; returns number densities
    END
