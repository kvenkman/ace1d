PRO ace_1d_major_solv, zmajnow, zmaj, pconst, model, mass, o_o2_terms, count
; Solve for major species 
; Refer to Dickinson, Ridley & Roble, 1984 / Sutton et al. 2015 for solver
; Setting up
; Source terms for the solver 
    s11 = o_o2_terms.s11
    s12 = o_o2_terms.s12
    s21 = o_o2_terms.s21
    s22 = o_o2_terms.s22
; Explicit source terms
    s10 = o_o2_terms.s10
    s20 = o_o2_terms.s20

    p00 = 1e6 ; 1e5 Pa, in dyne/cm^2
    p0  = model.p0 ; 5e-4, in dyne/cm^2
    t00 = 273.

; Refer to Dickinson, Ridley and Roble [1984] for these values
; These coefficients have a temperature & pressure dependence, but
; get cancelled out with the same dependence that D has in the calculation of phi_ij
    d0  = 0.2 ; in units of cm^2 s^-1    
;    d12 = 0.26
;    d13 = 0.18
;    d23 = 0.26
;    d21 = d12

    dt = deriv(model.zp, zmajnow.tn)
    dbarm = deriv(model.zp, zmajnow.barm)

    ; Derived parameters
    h0 = pconst.boltz*t00/(mass.n2*pconst.grav/pconst.avo) ; characteristic scale height of N2
    ;tau = exp(-model.zp)*(1.+fltarr(model.nlev))*p0*(h0^2.)/(p00*d0) ; diffusion time scale
;    h0 = pconst.boltz*zmajnow.tn/(mass.n2*pconst.grav/pconst.avo)
    tau = (1. + fltarr(model.nlev))*p0*(h0^2.)/(p00*d0) ; diffusion time scale ; 1.86e3; 
    ; d = d0*(p00/(model.p0*exp(-model.zp)))*(zmajnow.tn/t00)^1.75
; if count eq 10e3 then STOP
    sc_h = pconst.boltz*zmajnow.tn/(pconst.grav*zmajnow.barm/pconst.avo) ; Mean scale height in cm
    k_e = 5e-6*exp(-7-model.zp) ; Roble, 1987

	pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn) ; number density
	rho = zmajnow.barm*pkt/pconst.avo ; mass density in gm cm^-3
	drho = deriv(model.zp, rho)

    ; Calculating the alpha terms
    phi12 = 1.35; (d/d12)*(mass.n2/mass.o)
    phi13 = 1.11; (d/d13)*(mass.n2/mass.n2)
    phi21 = 0.673; (d/d21)*(mass.n2/mass.o2)
    phi23 = 0.769; (d/d23)*(mass.n2/mass.n2)

    alpha11 = (phi13 + (phi12 - phi13)*zmajnow.o) * (-1.)
    alpha22 = (phi23 + (phi21 - phi23)*zmajnow.o2) * (-1.)
    alpha12 = (phi12 - phi13)*zmajnow.o2 * (1.)
    alpha21 = (phi21 - phi23)*zmajnow.o * (1.)
    det_a = alpha11*alpha22 - alpha12*alpha21

    betam  = fltarr(2, 2, model.nlev)
    dbetam = fltarr(2, 2, model.nlev) ; Matrix containing derivative of elements of beta
    smatrix = fltarr(2, 2, model.nlev) ; Source term matrix
    smatrix_exp = fltarr(1, 2, model.nlev) ; Explicit source terms

    ; betam = inv(alpha)
    term1 = -((t00/zmajnow.tn)^0.25)*(zmajnow.barm/(tau*mass.n2))
    term1 = term1/det_a

	for i = 0, model.nlev - 1 do begin
		betam[*, *, i] = term1[i]*[[alpha22[i],-alpha12[i]], [-alpha21[i], alpha11[i]]]
		smatrix[*, *, i] = [[s11[i], s12[i]], [s21[i], s22[i]]]
		smatrix_exp[*, *, i] = [[s10[i]], [s20[i]]]
    endfor

    dbetam[0, 0, *] = deriv(model.zp, reform(betam[0, 0, *]))
    dbetam[1, 0, *] = deriv(model.zp, reform(betam[1, 0, *]))
    dbetam[0, 1, *] = deriv(model.zp, reform(betam[0, 1, *]))
    dbetam[1, 1, *] = deriv(model.zp, reform(betam[1, 1, *]))
    
    term2_1 = -(1. - mass.o2/zmajnow.barm - (dbarm/zmajnow.barm)) ; For O2
    term2_2 = -(1. -  mass.o/zmajnow.barm - (dbarm/zmajnow.barm)) ; For O

    term2  = fltarr(2, 2, model.nlev)
    term2[0, 0, *] = term2_1
    term2[1, 1, *] = term2_2

    dbeta_t2 = fltarr(2, 2, model.nlev)

;    for i = 0, model.nlev - 1 do dbeta_t2[*, *, i] = [[dtemp211[i], dtemp212[i]], [dtemp221[i], dtemp222[i]]]

    dbeta_t2[0, 0, *] = deriv(model.zp, term2_1*betam[0, 0, *])
    dbeta_t2[1, 0, *] = deriv(model.zp, term2_2*betam[1, 0, *])
    dbeta_t2[0, 1, *] = deriv(model.zp, term2_1*betam[0, 1, *])
    dbeta_t2[1, 1, *] = deriv(model.zp, term2_2*betam[1, 1, *])

    term3 = exp(-model.zp)*k_e
    term4 = dbarm/zmajnow.barm
    term5 = term3*term4
    term6 = exp(-model.zp)/model.del_time

    dterm3 = deriv(model.zp, term3)
    dterm5 = deriv(model.zp, term5)

    p     = fltarr(2,2,model.nlev)
    q     = fltarr(2,2,model.nlev)
    r     = fltarr(2,2,model.nlev)
    rhs   = fltarr(1,2,model.nlev)
; STOP
;;; The major species continuity equation is expressed as a block tridiagonal system of equations as O/O2 mmr are coupled
    FOR i = 0, model.nlev-1 DO BEGIN

        p[*, *, i] = $
                     betam[*, *, i]/(model.dz^2.) - $
                     (dbetam[*, *, i] + betam[*, *, i]##term2[*, *, i])/(2.*model.dz) + $
                     ; Eddy diffusion terms
                     term3[i]/(model.dz^2.)*[[1.,0.],[0.,1.]] - $
                     (dterm3[i] + term5[i])/(2.*model.dz)*[[1.,0.],[0.,1.]]

        q[*, *, i] = $
                     -2.*betam[*, *, i]/(model.dz^2.) + $
                     dbeta_t2[*, *, i] - $
                     ; (dbetam[*, *, i]##term2[*, *, i] + betam[*, *, i]##dterm2[*, *, i]) - $
                     ; Eddy diffusion terms
                     2.*term3[i]*[[1.,0.],[0.,1.]]/(model.dz^2.) + $
                     dterm5[i]*[[1.,0.],[0.,1.]] - $
                     ; Time derivative term
                     term6[i]*[[1.,0.],[0.,1.]] + $
                     ; Source term
                     exp(-model.zp[i])*smatrix[*, *, i]
                     
        r[*, *, i] = $
                     betam[*, *, i]/(model.dz^2.) + $
                     (dbetam[*, *, i] + betam[*, *, i]##term2[*, *, i])/(2.*model.dz) + $
                     ; Eddy diffusion terms
                     term3[i]/(model.dz^2.)*[[1.,0.],[0.,1.]] + $
                     (dterm3[i] + term5[i])/(2.*model.dz)*[[1.,0.],[0.,1.]]

        rhs[*, *, i] = $
                       ; Time derivative term
                       -term6[i]*[[zmajnow.o2[i]],[zmajnow.o[i]]] - $
                       ; Source term
                       exp(-model.zp[i])*smatrix_exp[*, *, i];[[s10[i]],[s20[i]]]

    ENDFOR
  
; Lower Boundary conditions

; This old set of boundary conditions define d[O]/dz = 0, and psi_o + psi_o2 = 0.234
; However, there are numerical issues with trying to define a local maxima of O density at the LB

;    m = p[0, 0, 0]/p[0, 1, 0]

;    q_o2_term = q[0, 0, 0] - m*q[0, 1, 0]
;    r_o2_term = r[0, 0 ,0] - m*r[0, 1, 0]

;    q_o_term = (q[1, 0, 0] - m*q[1, 1, 0]) - (2.*model.dz)*(p[1, 0, 0] - m*p[1, 1, 0]);*drho[0]/rho[0] ; *zmajnow.o[0]
;    r_o_term = (r[1, 0, 0] - m*r[1, 1, 0]) + (p[1, 0, 0] - m*p[1, 1, 0])
    
;    q[*, 0, 0] = [q_o2_term, q_o_term]; [0., 1.]; 
;    q[*, 1, 0] = [1., 1.]
    
;    r[*, 0, 0] = [r_o2_term, r_o_term]; [0., 0.]; 
;    r[*, 1, 0] = [0., 0.]
    
;    rhs[0, 0, 0] = rhs[0, 0, 0] - m*rhs[0, 1, 0]; 0.00754; 
;    rhs[0, 1, 0] = 0.234

;    p[*, *, 0] = fltarr(2, 2)
    
; Lower Boundary conditions
; This new set of LB simplifies the problem by assuming a F10.7 dependent constant 
; mass mixing ratio for O and O2

     q[*, 0, 0] = [0., 1.]
     q[*, 1, 0] = [1., 0.]
;    
     r[*, 0, 0] = [0., 0.]
     r[*, 1, 0] = [0., 0.]
    
     rhs[0, 0, 0] = zmajnow.o[0]
     rhs[0, 1, 0] = zmajnow.o2[0];0.234

     p[*, *, 0] = fltarr(2, 2)


;; Upper Boundary conditions, diffusive equilibrium for both O and O2
    p[*, *, model.nlev-1] = p[*, *, model.nlev-1] + r[*, *, model.nlev-1]
    ; There is a minus sign attached to the (2*model.dz) term below because a minus sign was previously absorbed into the expression for term2_i
    q[*, *, model.nlev-1] = q[*, *, model.nlev-1] - (2.*model.dz)*r[*, *, model.nlev-1]##[[term2_1[model.nlev-1], 0.], [0.,term2_2[model.nlev-1]]]
    r[*, *, model.nlev-1] = fltarr(2, 2)

;;; Block Tridiagonal matrix solver
    h     = fltarr(2,2,model.nlev)
    g     = fltarr(1,2,model.nlev)
    x     = fltarr(1,2,model.nlev)

;; Forward elimination
    h[*, *, 0] = mat_inv(q[*, *, 0])##r[*, *, 0]
    g[*, *, 0] = mat_inv(q[*, *, 0])##rhs[*, *, 0]
    
    FOR i = 1, model.nlev-1 DO BEGIN
        h[*, *, i] = mat_inv(q[*, *, i]-p[*, *, i]##h[*, *, i-1])##r[*, *, i]
        g[*, *, i] = mat_inv(q[*, *, i]-p[*, *, i]##h[*, *, i-1])##(rhs[*, *, i] - p[*, *, i]##g[*, *, i-1])
    ENDFOR

    o2_update = fltarr(model.nlev)
    o_update  = fltarr(model.nlev)

;; Backward substitution
    x[*, *, model.nlev-1] = g[*, *, model.nlev-1]
    o2_update[model.nlev-1] = x[0, 0, model.nlev-1] > 0.
    o_update[model.nlev-1] = x[0, 1, model.nlev-1] > 0.

    FOR i = model.nlev-2, 0, -1 DO BEGIN
        x[*, *, i] = (g[*, *, i] - h[*, *, i]##x[*, *, i+1])
        o2_update[i] = x[0, 0, i] > 0.
        o_update[i] = x[0, 1, i] > 0.
    ENDFOR

;;; N2 mass mixing ratios
    n2_update = (1. - o_update - o2_update) > 0.
;;; Updating mean mass, number densities
    barm_update = 1./(o_update/mass.o + o2_update/mass.o2 + n2_update/mass.n2)
    
    pkt = model.p0*exp(-model.zp)/(zmaj.tn*pconst.boltz)

    o_update_cm3  =  o_update*barm_update*pkt/mass.o
    o2_update_cm3 = o2_update*barm_update*pkt/mass.o2
    n2_update_cm3 = n2_update*barm_update*pkt/mass.n2

    zmaj.o = o_update
    zmaj.oden = o_update_cm3
    zmaj.o2 = o2_update
    zmaj.o2den = o2_update_cm3
    zmaj.n2 = n2_update
    zmaj.n2den = n2_update_cm3
    zmaj.barm = barm_update
    

END
