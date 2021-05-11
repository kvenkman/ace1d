FUNCTION ace_1d_oplus, x, prod, loss, zmajnow, zionnow, model, pconst
    ; Solving for O+ densities
    ; O+ is solved for using the equations for major ion diffusion from Schunk and Nagy [2001]
    ; eq. 5.54-5.59

    x_mass = 16. ; This solver only does O+
    tp = (zionnow.te + zionnow.ti)/2. ; Plasma temperature
    tr = (zmajnow.tn + zionnow.ti)/2. ; Reduced temperature
    theta_i = 75.*(!pi/180.); magnetic dip angle. Smithtro uses a value of 75 degrees in GAIT, while Roble uses a value of 90 degrees in the Global Mean Model.
    sin_sq_i = sin(theta_i)*sin(theta_i)
    c_in = [6.82, 6.64]*1e-10 ; Collision frequency coefficients, from Schunk and Nagy[2000], pg 105/107
    v_in_n2 = zmajnow.n2den*c_in[0]
    v_in_o2 = zmajnow.o2den*c_in[1]
    v_in_o  = zmajnow.oden*3.67e-11*sqrt(tr)*(1 - 0.064*alog10(tr))^2. > 0.

    v_in = (v_in_n2 + v_in_o2 + v_in_o); Assuming from Schunk [1988]; Eq. 30 

    pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn) ; number density
    h = pconst.boltz*zmajnow.tn/(pconst.grav*zmajnow.barm/pconst.avo) ; Mean scale height in cm
    hp = 2.*pconst.boltz*tp/(x_mass*pconst.grav/pconst.avo)
    d_a = 2.*pconst.boltz*tp/(x_mass*v_in/pconst.avo); ambipolar diffusion coefficient, Schunk and Nagy p.127
   
    lbc = (x[0] + prod[0]*model.del_time)/(1.+loss[0]*model.del_time)
    
    dtp = deriv(model.zp, tp)
    t1  = sin_sq_i*d_a/h
    t2 = (sin_sq_i*d_a/h)*(h/hp + dtp/tp)
    
    dt1 = deriv(model.zp, t1)
    dt2 = deriv(model.zp, t2)
    
    t3 = dt1 + t2
    
    ; Tridiagonal matrix elements    
    p   = fltarr(model.nlev)
    q   = fltarr(model.nlev)
    r   = fltarr(model.nlev)
    rhs = fltarr(model.nlev)

    p = (t1/(model.dz^2) - t3/(2.*model.dz))
    q = -2.*t1/(model.dz^2) +  dt2 - loss*h - h/model.del_time
    r = (t1/(model.dz^2) + t3/(2.*model.dz))
    rhs = (-h)*(x/model.del_time + prod)

    ; Lower boundary conditions (PCE)
    p[0] = 0.
    q[0] = 1.
    r[0] = 0.
    rhs[0] = lbc

    ; Upper boundary conditions (Diffusive Equilibrium)
    p[model.nlev-1] = p[model.nlev - 1] + r[model.nlev - 1]
    q[model.nlev-1] = q[model.nlev-1] - 2.*(model.dz)*(t2[model.nlev-1]/t1[model.nlev-1])*r[model.nlev-1]
    r[model.nlev-1] = 0.

    x_np = trisol(p, q, r, rhs) > 0.

    return, x_np
END
