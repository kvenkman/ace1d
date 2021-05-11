PRO ace_1d_dt, zmaj, zmajnow, zionnow, model, pconst, mass, heatterms, coolterms, cooltermse, cooltermsi

    sc_h = (pconst.boltz*zmajnow.tn*pconst.avo)/(zmajnow.barm*pconst.grav) ; Scale height
    nden = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn) ; number density
    rho = zmajnow.barm*nden/pconst.avo ; mass density
	
    k_e = 5e-6*exp(-7.-model.zp) ; Roble, 1987
    
    t1 = (zmajnow.kt/sc_h) + k_e*sc_h*zmajnow.cp*rho
    t2 = k_e*(sc_h^2.)*rho*pconst.grav
    
    t3 = model.p0/(pconst.grav*exp(model.zp))
    t4 = model.p0*zmajnow.cp/(pconst.grav*exp(model.zp))
    
    dt1 = deriv(model.zp, t1)
    dt2 = deriv(model.zp, t2)

    cool_implicit = (coolterms.no_cool_i + coolterms.co2_cool_i + coolterms.o3p_cool_i)
    cool_explicit = (coolterms.no_cool_e + coolterms.co2_cool_e + coolterms.o3p_cool_e + 0.*coolterms.nov)
    
    edep = heatterms.q_total + (cooltermse.neutrals_implicit*zionnow.te + cooltermsi.neut*zionnow.ti)
    edep_implicit = (cooltermse.neutrals_implicit + cooltermsi.neut)
  
    p = fltarr(model.nlev)
    q = fltarr(model.nlev)
    r = fltarr(model.nlev)
    rhs = fltarr(model.nlev)

    ; Coefficient matrix
    p = (t1/model.dz^2. - dt1/(2.*model.dz))
    q = (-1.)*(2.*t1/model.dz^2. + t4/model.del_time + t3*(cool_implicit + edep_implicit)/rho)
    r = (t1/model.dz^2. + dt1/(2.*model.dz))
    rhs = (-1.)*(t4*zmajnow.tn/model.del_time + (edep - cool_explicit)*t3/rho + dt2)

    ; Lower boundary condition, T = 177 K
    p[0] = 0;
    q[0] = 1;
    r[0] = 0;
    rhs[0] = zmajnow.tn[0];  177.
    
    ; Upper boundary condition, dT/dz = 0
    p[model.nlev-1] = p[model.nlev-1] + r[model.nlev-1]
    r[model.nlev-1] = 0.

    tn_update = trisol(p,q,r,rhs)
    zmaj.tn = tn_update
END
