FUNCTION ace_1d_nplus, x, prod, loss, zmajnow, zionnow, model, pconst
    ; Solving for N+ densities
    ; N+ is solved for using the equations for minor ion diffusion from Schunk and Nagy [2001]
    ; eq. 5.70-5.71

    x_mass = 14. ; This solver does N+
    tp = (zionnow.te + zionnow.ti)/2. ; Plasma temperature
    theta_i = 75.*(!pi/180.); magnetic dip angle. Smithtro uses a value of 75 degrees in GAIT, while Roble uses a value of 90 degrees in the Global Mean Model.
    sin_sq_i = sin(theta_i)*sin(theta_i)
    c_in = [7.47, 7.25, 4.42]*1e-10 ; Collision frequency coefficients, from Schunk and Nagy
    v_in_n2 = zmajnow.n2den*c_in[0]
    v_in_o2 = zmajnow.o2den*c_in[1]
    v_in_o  = zmajnow.oden*c_in[2]

    v_in = (v_in_n2 + v_in_o2 + v_in_o)

    pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn) ; number density
    h = pconst.boltz*zmajnow.tn/(pconst.grav*zmajnow.barm/pconst.avo) ; Mean scale height in cm
    hi = pconst.boltz*zionnow.ti/(x_mass*pconst.grav/pconst.avo)
    d_a = pconst.boltz*zionnow.ti/(x_mass*v_in/pconst.avo); minor ion diffusion coefficient, eq 5.71
   
    lbc = (x[0] + prod[0]*model.del_time)/(1.+loss[0]*model.del_time)
    
    dtp = deriv(model.zp, tp)
    t1  = sin_sq_i*d_a/h
    t2 = (sin_sq_i*d_a/h)*(h/hi + 2.*dtp/zionnow.ti) ; The 4th and 5th term of equation 5.72 are ignored here.

    dt1 = deriv(model.zp, t1)
    dt2 = deriv(model.zp, t2)
    
    t3 = dt1 + t2
    
    ; Tridiagonal matrix elements    
    p   = fltarr(model.nlev)
    q   = fltarr(model.nlev)
    r   = fltarr(model.nlev)
    rhs = fltarr(model.nlev)

    p = (t1/model.dz^2. - t3/(2.*model.dz))
    q = -2.*t1/(model.dz^2) +  dt2 - loss*h - h/model.del_time
    r = (t1/model.dz^2. + t3/(2.*model.dz))
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

    x_np = trisol(p, q, r, rhs)
    return, x_np
END
