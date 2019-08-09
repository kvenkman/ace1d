PRO ace_1d_km, zmaj, zminor, mass, model, pconst, k_m, flux

; Calculating molecular diffusion coefficients 

; This calculation for the diffusion coefficients was obtained from
; "Gaseous diffusion coefficients", Marrero and Mason [1972]
; Eq 2.2-8
; Also refer to Fundamentals of momentum, heat and mass transfer 5th ed. , (Eq 24-33) [Welty, 2000]
; Collision integrals from wright, bose, palmer, levin (2005)
; interpolated between temperatures given and temperatures in the model

;; N Collision integrals
    om_n_o=[8.32, 7.34,6.22, 5.26];  The collision integrals need to be in units of angstrom^2 (so DON'T convert to cm^-2)
    temp_n_o=[300.,500.,1000.,2000.]
    om_n_o_alt=interpol(om_n_o,temp_n_o,zmaj.tn)

    om_n_n2=[10.10, 8.57,7.70,6.65]
    temp_n_n2=[300.,600.,1000.,2000.]
    om_n_n2_alt=interpol(om_n_n2,temp_n_n2,zmaj.tn)

    om_n_o2=[7.56, 7.26, 6.55, 5.60]
    temp_n_o2=[500.,600.,1000.,2000.]
    om_n_o2_alt=interpol(om_n_o2,temp_n_o2,zmaj.tn)

;; NO Collision integrals
    om_no_o2=[11.39, 10.10, 9.75, 8.89, 7.74]
    temp_no_o2=[300.,500.,600.,1000.,2000.]
    om_no_o2_alt=interpol(om_no_o2,temp_no_o2,zmaj.tn)

    om_no_n2=[11.88, 10.61, 10.24, 9.35, 8.12]
    temp_no_n2=[300.,500.,600.,1000.,2000.]
    om_no_n2_alt=interpol(om_no_n2,temp_no_n2,zmaj.tn)

    om_no_o=[7.57, 7.27, 6.55, 5.62]
    temp_no_o=[500.,600.,1000.,2000.]
    om_no_o_alt=interpol(om_no_o,temp_no_o,zmaj.tn)

;; O Collision integrals
    om_o_o2 = [9.10, 7.58, 6.74, 5.7]
    temp_o_o2 = [300., 600., 1000., 2000.]
    om_o_o2_alt = interpol(om_o_o2, temp_o_o2, zmaj.tn)

    om_o_n2 = [8.07, 5.93, 5.17]
    temp_o_n2 = [300., 1000., 2000.]
    om_o_n2_alt = interpol(om_o_n2, temp_o_n2, zmaj.tn)

;; O2 Collision integrals
    om_o2_o = om_o_o2
    temp_o2_o = temp_o_o2
    om_o2_o_alt = om_o_o2_alt

    om_o2_n2 = [10.16, 7.39, 6.42]
    temp_o2_n2 = [300., 1000., 2000.]
    om_o2_n2_alt = interpol(om_o2_n2, temp_o2_n2, zmaj.tn)

;; N2 Collision integrals
    om_n2_o = om_o_n2
    temp_n2_o = temp_o_n2
    om_n2_o_alt = om_o_n2_alt

    om_n2_o2 = om_o2_n2
    temp_n2_o2 = temp_o2_n2
    om_n2_o2_alt = om_o2_n2_alt

;; Reduced masses
    mu_no = mass.no*[mass.n2/(mass.n2+mass.no), mass.o2/(mass.o2+mass.no), mass.o/(mass.o+mass.no)]
    mu_n = mass.n4s*[mass.n2/(mass.n2+mass.n4s), mass.o2/(mass.o2+mass.n4s), mass.o/(mass.o+mass.n4s)]
    mu_o = mass.o*[mass.n2/(mass.n2+mass.o), mass.o2/(mass.o2+mass.o), mass.o/(mass.o+mass.o)]
    mu_o2 = mass.o2*[mass.n2/(mass.n2+mass.o2), mass.o2/(mass.o2+mass.o2), mass.o/(mass.o+mass.o2)]
    mu_n2 = mass.n2*[mass.n2/(mass.n2+mass.n2), mass.o2/(mass.o2+mass.n2), mass.o/(mass.o+mass.n2)]

    oden  = zmaj.oden
    o2den = zmaj.o2den
    n2den = zmaj.n2den

    pkt = model.p0*exp(-model.zp)/(zmaj.tn*pconst.boltz); Total density
    p_atm = model.p0*exp(-model.zp)/1.01325e6 ; Converting pressure from dyne/cm^2 into atmospheres
    rho = zmaj.barm*pkt/pconst.avo
    h = pconst.boltz*zmaj.tn/(pconst.grav*zmaj.barm/pconst.avo) ; Mean scale height in cm
    k_e = 5e-6*exp(-7-model.zp)*(h^2.) ; Roble, 1987
       
    dt = deriv(model.zp, zmaj.tn)/h
    
;;; D_NO,j
    d_no_o=(1./!pi)*0.008258*sqrt(1./(2.*mu_no[2]))*(zmaj.tn^1.5)/(p_atm*om_no_o_alt)
    d_no_o2=(1./!pi)*0.008258*sqrt(1./(2.*mu_no[1]))*(zmaj.tn^1.5)/(p_atm*om_no_o2_alt)
    d_no_n2=(1./!pi)*0.008258*sqrt(1./(2.*mu_no[0]))*(zmaj.tn^1.5)/(p_atm*om_no_n2_alt)

;    d_no_o  = 4.5e17*sqrt(zmaj.tn)/pkt ; From Strobel, 1971
;    d_no_o2 = 2.9e17*sqrt(zmaj.tn)/pkt
;    d_no_n2 = 2.9e17*sqrt(zmaj.tn)/pkt
    
;;; D_N,j
    d_n_o=(1./!pi)*0.008258*sqrt(1./(2.*mu_n[2]))*(zmaj.tn^1.5)/(p_atm*om_n_o_alt)
    d_n_o2=(1./!pi)*0.008258*sqrt(1./(2.*mu_n[1]))*(zmaj.tn^1.5)/(p_atm*om_n_o2_alt)
    d_n_n2=(1./!pi)*0.008258*sqrt(1./(2.*mu_n[0]))*(zmaj.tn^1.5)/(p_atm*om_n_n2_alt)

;    d_n_o  = 6.9e17*sqrt(zmaj.tn)/pkt ; From Strobel, 1971
;    d_n_o2 = 4.6e17*sqrt(zmaj.tn)/pkt
;    d_n_n2 = 4.5e17*sqrt(zmaj.tn)/pkt

;;; D_O,j
    ;d_o_o=(1./!pi)*0.008258*sqrt(1./(2.*mu_o[2]))*(zmaj.tn^1.5)/(p_atm*om_o_o_alt)
    d_o_o2=(1./!pi)*0.008258*sqrt(1./(2.*mu_o[1]))*(zmaj.tn^1.5)/(p_atm*om_o_o2_alt)
    d_o_n2=(1./!pi)*0.008258*sqrt(1./(2.*mu_o[0]))*(zmaj.tn^1.5)/(p_atm*om_o_n2_alt)

;;; D_O2,j
    d_o2_o=(1./!pi)*0.008258*sqrt(1./(2.*mu_o2[2]))*(zmaj.tn^1.5)/(p_atm*om_o2_o_alt)
    ;d_o2_o2=(1./!pi)*0.008258*sqrt(1./(2.*mu_o2[1]))*(zmaj.tn^1.5)/(p_atm*om_o2_o2_alt)
    d_o2_n2=(1./!pi)*0.008258*sqrt(1./(2.*mu_o2[0]))*(zmaj.tn^1.5)/(p_atm*om_o2_n2_alt)

;;; D_N2,j
    d_n2_o=(1./!pi)*0.008258*sqrt(1./(2.*mu_n2[2]))*(zmaj.tn^1.5)/(p_atm*om_n2_o_alt)
    d_n2_o2=(1./!pi)*0.008258*sqrt(1./(2.*mu_n2[1]))*(zmaj.tn^1.5)/(p_atm*om_n2_o2_alt)
    ;d_n2_n2=(1./!pi)*0.008258*sqrt(1./(2.*mu_n2[0]))*(zmaj.tn^1.5)/(p_atm*om_n2_n2_alt)
    
    k_m.n4s = pkt*(n2den/d_n_n2 + o2den/d_n_o2 + oden/d_n_o)^(-1.)
                 
    k_m.no = pkt*(n2den/d_no_n2 + o2den/d_no_o2 + oden/d_no_o)^(-1.)
          
    k_m.o = pkt*(n2den/d_o_n2 + o2den/d_o_o2)^(-1.)
    
    k_m.o2 = pkt*(n2den/d_o2_n2 + oden/d_o2_o)^(-1.)
                
    k_m.n2 = pkt*(o2den/d_n2_o2 + oden/d_n2_o)^(-1.)
    
;    IF(count GT 1000) THEN BEGIN
                    
    dn_o  = deriv(model.zp, zmaj.oden)
    dn_o2 = deriv(model.zp, zmaj.o2den)
    dn_n2 = deriv(model.zp, zmaj.n2den)

    d_o = (zmaj.o2den/(pkt*d_o_o2)  + zmaj.n2den/(pkt*d_o_n2))^(-1.)
    d_o2 = (zmaj.oden/(pkt*d_o2_o)  + zmaj.n2den/(pkt*d_o2_n2))^(-1.)
    d_n2 = (zmaj.oden/(pkt*d_n2_o)  + zmaj.o2den/(pkt*d_n2_o2))^(-1.)

    phi_o =  -d_o*( $
              dn_o/h + zmaj.oden*mass.o*pconst.grav/(pconst.boltz*zmaj.tn*pconst.avo) + $
              zmaj.oden*dt/zmaj.tn - zmaj.oden*(flux.o2/(pkt*d_o_o2) + flux.n2/(pkt*d_o_n2))) - $
              k_e*(dn_o/h + zmaj.oden/h + zmaj.oden*dt/zmaj.tn)
   
    phi_o2 =  -d_o2*( $
              dn_o2/h + zmaj.o2den*mass.o2*pconst.grav/(pconst.boltz*zmaj.tn*pconst.avo) + $
              zmaj.o2den*dt/zmaj.tn - zmaj.o2den*(flux.o/(pkt*d_o2_o) + flux.n2/(pkt*d_o2_n2))) - $
              k_e*(dn_o2/h + zmaj.o2den/h + zmaj.o2den*dt/zmaj.tn)
   
    phi_n2 =  -d_n2*( $
              dn_n2/h + zmaj.n2den*mass.n2*pconst.grav/(pconst.boltz*zmaj.tn*pconst.avo) + $
              zmaj.n2den*dt/zmaj.tn - zmaj.n2den*(flux.o/(pkt*d_n2_o) + flux.o2/(pkt*d_n2_o2))) - $
              k_e*(dn_n2/h + zmaj.n2den/h + zmaj.n2den*dt/zmaj.tn)
    
    
    ; Not the fluxes of NO and N(4S), rather the frictional terms that appear in their diffusion expressions.
    phi_no  = flux.o/(pkt*d_no_o) + flux.o2/(pkt*d_no_o2) + flux.n2/(pkt*d_no_n2)
    phi_n4s = flux.o/(pkt*d_n_o) + flux.o2/(pkt*d_n_o2) + flux.n2/(pkt*d_n_n2)
    
    flux.o = phi_o
    flux.o2 = phi_o2
    flux.n2 = phi_n2
    flux.no = phi_no
    flux.n4s = phi_n4s
    
;    ENDIF
END
