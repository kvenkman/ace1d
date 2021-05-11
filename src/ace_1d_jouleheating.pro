PRO ace_1d_jouleheating, zmajnow, zminornow, zionnow, pconst, mass, model, heatterms, model_sun, count

    ; Roble '87 used 5.7mV/m - 3.6mV/m as the range of E-field
    ; values for Joule heating for F10.7 between 67 and 243

    ; Calculations from this model agree with the above values for the E-field
    
    ; We assume the hemispheric energy input due to Joule heating is 95 GW
    ; From Knipp [2004]. Foster [1983], and Roble [1987] assumed 70 GW
    
    ; For a variable Joule Heating rate
    ; 16e10 for f10.7 = 67
    ; 10e10 for f10. = 250
    ;
    ;m = (10e10-16e10)/(250.-67.)
    ;c =  10e10 - m*250.
    ; 
    ;h_input = m*model_sun.f107d + c ; Hemispheric energy input, Watts
    ;;;;;;
    
    h_input = 7e10 ; Hemispheric energy input, Watts
    
    h_area = 4.*!pi*(6371e3)^2. ; Area in m^2
    q_joule_flux = h_input/h_area ; height integrated energy flux, W/m^2
    
    dz = model.dz*(pconst.boltz*zmajnow.tn)/(pconst.grav*zmajnow.barm/pconst.avo)*1e-2 ; meters
    
	b_tesla = 0.4e-4 ; B-field in Tesla

    tr = (zmajnow.tn + zionnow.ti)/2.
	
	; Collision frequencies, Schunk & Nagy. Pg. 105-107
    v_op_o   = 3.67e-11*zmajnow.oden*(tr^0.5)*(1.-0.064*alog10(tr))^2
    v_op_o2  = 6.64e-10*zmajnow.o2den
    v_op_n2  = 6.82e-10*zmajnow.n2den
    v_o2p_o  = 2.31e-10*zmajnow.oden
    v_o2p_o2 = 2.59e-11*zmajnow.o2den*(tr^0.5)*(1.-0.073*alog10(tr))^2
    v_o2p_n2 = 4.13e-10*zmajnow.n2den
    v_nop_o  = 2.44e-10*zmajnow.oden
    v_nop_o2 = 4.27e-10*zmajnow.o2den
    v_nop_n2 = 4.34e-10*zmajnow.n2den
    
    v_o2p = v_o2p_o2 + v_o2p_o + v_o2p_n2
    v_op  = v_op_o2 + v_op_o + v_op_n2
    v_nop = v_nop_o2 + v_nop_o + v_nop_n2
    
    ; pg.109, pg. 141 (Eq. 5.114)
    v_en = 2.33e-11*zmajnow.n2den*zionnow.te*(1. - 1.21e-4*zionnow.te) + $
           1.82e-10*zmajnow.o2den*(zionnow.te^0.5)*(1. + 3.6e-2*(zionnow.te^0.5)) + $
           8.9e-11*zmajnow.oden*(zionnow.te^0.5)*(1. + 5.7e-4*zionnow.te)
    
    ; gyrofrequencies are given as qB/m
    w_o2p = pconst.q_e*b_tesla*pconst.avo/(1e-3*mass.o2)
    w_op = pconst.q_e*b_tesla*pconst.avo/(1e-3*mass.o)
    w_nop = pconst.q_e*b_tesla*pconst.avo/(1e-3*mass.no)
    w_e = pconst.q_e*b_tesla/(1e-3*pconst.m_e)
    
    r_o2p = v_o2p/w_o2p
    r_op = v_op/w_op
    r_nop = v_nop/w_nop
    r_e = v_en/w_e
    
    qfac = pconst.q_e/b_tesla
    ; Schunk and Nagy, pg. 141 (Eq 5.117)
    ; Pedersen conductivity, S/m
    s_ped = 1e6*qfac*($          ; 1e6 is for conversion of densities in cm^-3 -> m^-3
             zionnow.o_p*r_op/(1. + r_op*r_op)        + $
    		 zionnow.o2_p*r_o2p/(1. + r_o2p*r_o2p)    + $
    		 zionnow.no_p*r_nop/(1. + r_nop*r_nop)    + $
    		 zionnow.e*r_e/(1. + r_e*r_e))

    efield = sqrt(q_joule_flux/total(s_ped*dz))

    q_joule = s_ped*efield*efield*10.; * 1e-6 * 1e7 ; Conversion of J s^-1 m^-3 to ergs s^-1 cm^-3

    heatterms.s_ped = s_ped ; In S/m
	heatterms.q_joule = q_joule
	heatterms.q_total = heatterms.q_total + q_joule

END
