PRO ace_1d_cooling, zmajnow, zminornow, zionnow, termsno, pconst, model, coolterms
; This routine might eventually simply call more complex NO and CO2 cooling procedures

; Collisional cooling, form given by Kockarts [1980]
; NO Cooling
; There are many different rates coefficients that are used for the NO + O collisional excitation rate,
; starting from 6.5E-11 [Fernando & Smith] used by Kockarts et al [1980], 
; 4.2E-11 [Hwang et. al] used by Qian et al., and the most 
; recent value of 2.1E-11 calculated by Hwang et. al

; Hwang et al. value:
; k_01 = 4.2e-11
; Fernando & Smith value:
; k_01 = 6.5e-11
; Caridade et al value:
  k_01 = 2.1e-11

a_10 = 13.3 ; Einstein coefficient for NO (v = 1-0) emission
lambda_no = 5.3e-4
photon_53 = pconst.h*pconst.c/lambda_no ; Energy of 5.3 micron photon

; "Dilution" factor. Ref: Kockarts 1980.
omega = k_01*zmajnow.oden/(k_01*zmajnow.oden + a_10)

; The NO cooling rate due to collisional excitation in ergs cm^-3 s^-1
no_cool_a = photon_53*a_10*zminornow.no*omega
no_cool_b = photon_53/pconst.boltz

; Implicit and explicit terms are needed for the temperature solver
no_cool_explicit = no_cool_a*exp(-no_cool_b/zmajnow.tn)*(1.-no_cool_b/zmajnow.tn)
no_cool_implicit = no_cool_a*no_cool_b*exp(-no_cool_b/zmajnow.tn)/(zmajnow.tn^2.)
; But this is the expression for the NO cooling at the current time step
no_cool = no_cool_a*exp(-no_cool_b/zmajnow.tn)

; Calculating chemiluminescent NO emission
nov = ace_1d_nov(zmajnow, zminornow, zionnow, termsno, model, pconst)

;;;;;;;;;;;;;;;;;;;;;;;;
; CO2 cooling
; Dickinson [1984], Gustav M Shved [2003], Castle [2012]
a_co2 = fltarr(model.nlev)
b_co2 = fltarr(model.nlev)

lambda_co2 = 15e-4
photon_15 = pconst.h*pconst.c/lambda_co2
;a10_co2 = 

; a(co2) [Dickinson, 1984]
; [This is cooling by VT collisions with O2/N2]
mask1 = (zmajnow.tn LT 200.)
mask2 = (zmajnow.tn GE 200.)
a_co2 = 2.5e-15*mask1 + 2.5e-15*(1.+0.03*(zmajnow.tn-200.))*mask2

; b(co2) [Shved, 2003]
;mask1 = zmajnow.tn LT 260.  				    
;mask2 = (zmajnow.tn GE 260. AND zmajnow.tn LE 300.)
;mask3 = zmajnow.tn GT 300.
;b_co2 = mask1*1.56e-12 + mask2*(2.6-0.004*zmajnow.tn)*1.0e-12 + mask3*1.4e-12

; b(co2) [Castle, 2012]
; [This is cooling by VT collisions with O]
mask1 = (zmajnow.tn LE 490)
mask2 = (zmajnow.tn GT 490)

b_co2 = mask1*(3.8 - 9.51e-3*zmajnow.tn + 9.32e-6*(zmajnow.tn^2))*1e-12 + $
	    mask2*1.38e-12 ;Assuming constant k for T> 490 K

; CO2 cooling at current timestep
co2_cool = 2.65e-13*zminornow.co2*exp(-960./zmajnow.tn)*((zmajnow.o2den + zmajnow.n2den)*a_co2 + b_co2*zmajnow.oden)

; Implicit and explicit terms for temperature solver
co2_term1 = 2.65e-13*zminornow.co2*((zmajnow.o2den + zmajnow.n2den)*a_co2 + b_co2*zmajnow.oden)
co2_term2 = pconst.h*(pconst.c/lambda_co2)/pconst.boltz

co2_cool_explicit = co2_term1*exp(-co2_term2/zmajnow.tn)*(1.-co2_term2/zmajnow.tn)
co2_cool_implicit = co2_term1*exp(-co2_term2/zmajnow.tn)*co2_term2/(zmajnow.tn^2.)

;;;;;;;;;;;;;;;;;;;;
; O(3P) cooling
;
; The cooling is given by Eq 4. of Kockarts and Peetermans [1970], multplied by a 'masking factor' accounting 
; for the varying optical thickness of the atmosphere at 63 microns. However, this calculation should 
; be revisited and implemented more rigorously.

; From TIE-GCM newton.F
mfac =     [ $  
           0.1000E-01, 0.1000E-01, 0.1000E-01, 0.1000E-01, 0.1000E-01, $
           0.3000E-01, 0.5000E-01, 0.7500E-01, 0.1000E+00, 0.1500E+00, $
           0.2000E+00, 0.3000E+00, 0.4000E+00, 0.4750E+00, 0.5500E+00, $
           0.6250E+00, 0.7000E+00, 0.7250E+00, 0.7500E+00, 0.7750E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, 0.8000E+00, $
           0.8000E+00, 0.8000E+00]


lambda_o3p = 63e-4

; O(3P) cooling at the current timestep
o3p_cool = mfac*1.69e-18*zmajnow.oden*exp(-228./zmajnow.tn)/(1. + 0.6*exp(-228./zmajnow.tn) + 0.2*exp(-326./zmajnow.tn))

; The following comes out of a taylor series expansion of the O(3P) cooling expression about T(n)
hvk = pconst.h*(pconst.c/lambda_o3p)/pconst.boltz
o3p_a = exp(-hvk/zmajnow.tn)
o3p_b = 1. + 0.6*exp(-228./zmajnow.tn) + 0.2*exp(-326./zmajnow.tn)

; Implicit and explicit cooling for the solver
o3p_cool_implicit = 1.69e-18*zmajnow.oden*mfac*((o3p_a/o3p_b)*(228./zmajnow.tn^2.) - (o3p_a/o3p_b^2.)*(0.6*exp(-228./zmajnow.tn)*$
                    (228./zmajnow.tn^2.) + 0.2*exp(-326./zmajnow.tn)*(326./zmajnow.tn^2.)))
o3p_cool_explicit = 1.69e-18*zmajnow.oden*mfac*(o3p_a/o3p_b - zmajnow.tn*((o3p_a/o3p_b)*(228./zmajnow.tn^2.) - (o3p_a/o3p_b^2.)*$
                    (0.6*exp(-228./zmajnow.tn)*(228./zmajnow.tn^2.) + 0.2*exp(-326./zmajnow.tn)*(326./zmajnow.tn^2.))))


;o3p_t1 = zmajnow.oden * (pconst.h*(pconst.c/lambda_o3p)) * 8.95e-5 * 0.6 ; Last two factors are the Einstein coefficient and the statistical weight for the O(3P1) level

;;;;;;;;;;;;;;;;;;;;;;;;;;
;Bringing it all together

net_cool = o3p_cool + co2_cool + nov + no_cool

coolterms = {noc:no_cool, nov:nov, o3p_cool:o3p_cool, co2_cool:co2_cool, $
             no_cool_i:no_cool_implicit, no_cool_e:no_cool_explicit, $
             co2_cool_i:co2_cool_implicit, co2_cool_e:co2_cool_explicit, $
             o3p_cool_i:o3p_cool_implicit, o3p_cool_e:o3p_cool_explicit, $
             net_cool:net_cool}

END
