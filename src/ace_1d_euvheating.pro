PRO ace_1d_euvheating, zmajnow, zcol, xsec, branching, solspec, heatterms, pconst, model, model_sun

; Heating due to direct impact of fast photoelectrons
; The 5% efficiency factor is obtained from Roble [1987]
    lambda = 0.5*(solspec.wave1 + solspec.wave2)*1e-8 ; in centimeters
    euv_eff = .05 ; EUV heating efficiency due to photoelectrons
    
    euv_heating = fltarr(model.nlev)
    
    i0 = where(solspec.wave2 le 1050.0) ; only looking at EUV, not FUV
    FOR i = i0[0], solspec.n_wave - 1 DO BEGIN
        photon_e = pconst.h*pconst.c/lambda[i] ; in ergs
        
        ; Optical Depth
        tau = zcol.so  * xsec.abs(0, i) + $
              zcol.so2 * xsec.abs(1, i) + $
              zcol.sn2 * xsec.abs(2, i)
        ; Solar flux at each altitude at a given wavelength, in units of photons cm^-2 s^-1
        f = solspec.current[i]*exp(-tau)
        
        ; Absorption rate of photons at each altitude at a given wavelength, in units of photons cm^-3 s^-1 
        absorption = f*(xsec.abs(0, i)*zmajnow.oden  + $ 
                        xsec.abs(1, i)*zmajnow.o2den + $
                        xsec.abs(2, i)*zmajnow.n2den   $
                       )
                       
        euv_heating = euv_heating + (absorption*photon_e) ; converting # of photons into energy in ergs
 
    ENDFOR
   heatterms.q_euv = euv_heating*euv_eff ; Finally in units of ergs cm^-3 s^-1
   
; Schumann-Runge bands
    srb_heating = fltarr(model.nlev)
    srb_afac = 0.67     ; shumann-runge band heating parameters [Strobel, 1978]
    srb_bfac = 3.44e9    
    srb_cfac = 2.43e-19  
    
    loc1 = (zcol.so2 GE 1e18)
    loc2 = (zcol.so2 LT 1e18)

    srb_heating = loc1*(1./(srb_afac*zcol.so2 + srb_bfac*sqrt(zcol.so2))) + $
                  loc2*srb_cfac
    
    srb_heating = srb_heating*zmajnow.o2den
    
    heatterms.q_srb = srb_heating
    
; Schumann-Runge continuum
    src_heating = fltarr(model.nlev)
    
    do2 = 5.15 ; dissociation energy of O2
    e630 = pconst.h*pconst.c/(630e-7*pconst.ev2erg) ; Energy of 630 nm photon (relaxtion of O(1D) to ground)

    FOR i = 0, solspec.n_wave-1 DO BEGIN ; Loop over wavelength, but really concerned with between 1750 - 1050 A
        if solspec.wave1[i] ge 1200. then begin
            photon_e = pconst.h*pconst.c/lambda[i] ; in ergs
            ec = (do2+e630) - (1.0-branching.pdo2_o1dyield[i])*e630 ; Equation 5 of DeMajistre et al 2001
            src_heating = src_heating + ((photon_e - ec*1.6e-12)>0.0)*zmajnow.o2den*xsec.abs(1, i)*solspec.zflux[i,*]  ; 1.6e-12 converts eV into ergs
        endif
    ENDFOR
    
    heatterms.q_src = src_heating
    
; Total heating
    heatterms.q_total = heatterms.q_total + (heatterms.q_euv + heatterms.q_srb + heatterms.q_src)
    
END
