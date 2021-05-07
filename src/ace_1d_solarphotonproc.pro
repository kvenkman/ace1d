PRO ace_1d_solarphotonproc, zmaj, model, pconst, solspec, zcol, zpid, zei, zion, heattermse, edep, branching, xsec, model_sun, count

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;

;
;OUTPUT:	I_O-	Ionization freq for atomic oxygen.
;			I_O2-	Ionization freq for molecular oxygen
;			I_N2-	Ionization freq for molecular nitrogen
;			I_NO	Ionization freq for nitric oxide
;			I_N	Ionization freq for atomic nitrogen
;			J_N2-	Photolysis rate for molecular nitrogen
;			J_NO	Photolysis rate for nitric oxide
;			DI_N2 	Dissociative ionization rate for N2
;			DI_O2	Dissociative ionization rate for O2
;			PEI_N2	Photoelectron impact ionization of N2
;			PEI_O2	Photoelectron impact ionization of O2
;			PEI_O	Photoelectron impact ionization of O
;     PED_N2  Photeelectron impact excitation efficiency
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
;@ace_chem1d_common_blocks.pro
 
;****************** CALC OPTICAL DEPTH MATRICIES   ************************
tau_o  =  reform(xsec.abs[0,*])#zcol.so 
tau_o2 =  reform(xsec.abs[1,*])#zcol.so2
tau_n2 =  reform(xsec.abs[2,*])#zcol.sn2
tau    =  ((tau_o + tau_o2 + tau_n2) * (-1.0));>(-55.)
flux=solspec.current
; save the altitude dependent solar flux for use in later calculations (O(1D) production etc.)
for i=0,solspec.n_wave-1 do solspec.zflux[i,*]=flux[i]*exp(tau[i,*])

; br_pi = photoionization branching ratio
sigi_o =reform(xsec.abs[0,*]*xsec.br_pi[0,*]) ; Ionization cross sections
sigi_o2=reform(xsec.abs[1,*]*xsec.br_pi[1,*])
sigi_n2=reform(xsec.abs[2,*]*xsec.br_pi[2,*])

;************************ CALC IONIZATION FREQUENCIES **************************
zpid.i_o     =(sigi_o*                     flux)#exp(tau)
zpid.i_o_op2p=(sigi_o*branching.pio_op_2p *flux)#exp(tau)
zpid.i_o_op2d=(sigi_o*branching.pio_op_2d *flux)#exp(tau)
zpid.i_o_op4s=(sigi_o*branching.pio_op_4s *flux)#exp(tau)

;zpid.i_o2 = (sigi_o2*flux)#exp(tau)
;zpid.i_n2 = (sigi_n2*flux)#exp(tau)
zpid.i_o2 = (xsec.abs[1,*]*(reform(xsec.br_pi[1,*] - branching.pdio2))*flux)#exp(tau)
zpid.i_n2 = (xsec.abs[2,*]*(reform(xsec.br_pi[2,*] - branching.pdin2))*flux)#exp(tau)
zpid.i_n4s= (xsec.sigin4s*flux)#exp(tau)

; photodissociation or O2, we save the J_O2 by wavelength to use later in the calculation of O(1D) production
;zpid.j_o2 = (xsec.abs[1,*]*branching.pdo2*(2.-branching.pdo2_o1dyield)*flux)#exp(tau)

zpid.j_o2 = (xsec.abs[1,*]*branching.pdo2*flux)#exp(tau)

for i=0,n_elements(zpid.j_o2wvln[1,*])-1 do $ 
	zpid.j_o2wvln[*,i]=(xsec.abs[1,i]* branching.pdo2[i]*flux[i])*exp(tau)*zmaj.o2den

zpid.j_o2_o1d = (xsec.abs[1,*]*branching.pdo2*branching.pdo2_o1dyield*flux)#exp(tau)

; photodissociation of N2
zpid.j_n2 = (xsec.abs[2,*]* branching.pdn2 *flux)#exp(tau)

;	I_N20 = TOTAL(SIGI_N2*FLUX)

; the following needs to be improved; don't see any value for Lyman alpha absorption of NO in the arrays. 
lyman=solspec.current[11]	;solar ly-alpha
sigi_no=2.02e-18				;ly-a x-section for no
sig_ly_o2=xsec.abs[1,11]					;ly-a x-section for o2
zpid.i_no=(lyman*sigi_no*exp(-1.*sig_ly_o2*zcol.so2))
;0.*
;2.91E11*(1.+0.2*(model_sun.f107d-65.)/100.)*2.E-18*exp(-8.E-21*zcol.so2)
zpid.di_o2=(xsec.abs[1,*]*branching.pdio2*flux)#exp(tau)
zpid.di_n2=(xsec.abs[2,*]*branching.pdin2*flux)#exp(tau)

; Photoleectron

zei.pei_o  = (sigi_o  * reform(xsec.pepiscale[0,*]) * flux)#exp(tau)
zei.pei_o2 = (sigi_o2 * reform(xsec.pepiscale[1,*] - xsec.pepiscale_pedio2) * flux)#exp(tau)
zei.pei_n2 = (sigi_n2 * reform(xsec.pepiscale[2,*] - xsec.pepiscale_pedin2) * flux)#exp(tau)

;zei.pei_o_op2p = (sigi_o  * reform(xsec.pepiscale[0,*]* branching.peio_op2p) * flux)#exp(tau)
;zei.pei_o_op2d = (sigi_o  * reform(xsec.pepiscale[0,*]* branching.peio_op2d) * flux)#exp(tau)
;zei.pei_o_op4s = (sigi_o  * reform(xsec.pepiscale[0,*]* branching.peio_op4s) * flux)#exp(tau)

zei.pei_o_op2p = (sigi_o  * reform(branching.peio_op2p) * flux)#exp(tau)
zei.pei_o_op2d = (sigi_o  * reform(branching.peio_op2d) * flux)#exp(tau)
zei.pei_o_op4s = (sigi_o  * reform(branching.peio_op4s) * flux)#exp(tau)

zei.di_o2=(sigi_o2*xsec.pepiscale_pedio2*flux)#exp(tau)
zei.di_n2=(sigi_n2*xsec.pepiscale_pedin2*flux)#exp(tau)

zei.ped_o2 = (sigi_o2 *xsec.pepiscale_pedo2 *flux)#exp(tau)
zei.ped_n2 = (sigi_n2 *xsec.pepiscale_pedn2 *flux)#exp(tau)

; heating of ambient electrons by photoelectrons
; using the method of Smithro and Solomon 2008, full parameterization
nc=7
c_0_55   = [1.468, 9.229e-1,  4.956e-2, -1.897e-2, -3.934e-3, -2.643e-4, -5.980e-6]
c_55_105 = [1.020, 1.540e-2, -6.858e-3, -8.528e-3, -2.052e-3, -1.634e-4, -4.314e-6]
avee0_55 = [30.44, 40.74, 48.97, 68.00, 125.9, 216.0, 442.6, 1119., 2325., 4658.]

r=zion.e/(zmaj.n2den + zmaj.o2den + zmaj.oden) ;>1e-6

; need ionization rate from photons less than and greater than 55 nm
ind0_55  =where(solspec.wave2 le 550.) ; solar flux is in angstroms 
ind55_105=where(solspec.wave1 gt 550.)
mask0_55  =fltarr(solspec.n_wave)
mask55_105=fltarr(solspec.n_wave)
mask0_55[ind0_55]    =avee0_55 ; see Equation 7 of Smithro and Solomon
mask55_105[ind55_105]=1.0
flux0_55  =flux*mask0_55
flux55_105=flux*mask55_105

exptausbo=tau & exptausbo2=tau & exptausbn2=tau
for i=0,solspec.n_wave-1 do begin
  exptausbo[i,*] =exp(exptausbo[i,*]) *zmaj.oden
  exptausbo2[i,*]=exp(exptausbo2[i,*])*zmaj.o2den
  exptausbn2[i,*]=exp(exptausbn2[i,*])*zmaj.n2den
end

ionzr0_55       =(sigi_o*flux0_55)  #exptausbo + (sigi_o2*flux0_55)  #exptausbo2 + (sigi_n2*flux0_55)  #exptausbn2
ionzr55_105     =(sigi_o*flux55_105)#exptausbo + (sigi_o2*flux55_105)#exptausbo2 + (sigi_n2*flux55_105)#exptausbn2
; heat efficiency, epsilon
eps0_55   = r*0.0
eps55_105 = r*0.0
for i=0,nc-1 do begin
eps0_55   = eps0_55   +   c_0_55[i]*(alog(r)^i)
eps55_105 = eps55_105 + c_55_105[i]*(alog(r)^i)  
endfor
eps0_55  =exp(eps0_55)
eps55_105=exp(eps55_105)

; heating rate
heattermse.q_pe0_55   = eps0_55*ionzr0_55 * pconst.ev2erg
heattermse.q_pe55_105 = eps55_105*ionzr55_105 * pconst.ev2erg
;if count EQ 4998 then STOP
; Chemical production due to photon processes:

;O(1D) from dissociation of O2
; zpid.j_o2 = (xsec.abs[1,*]* branching.pdo2*branching.pdo2_o1dyield *flux)#exp(tau)
;coeff = 5.1E-11*(300./zmaj.tn)^(1.16) ; Rate coefficient
;ace_1d_update_crh_solar, spindex.n2_p, spindex.o2, zion.n2_p, zmaj.o2den, $
;                   coeff, exomatrix[spindex.n2_p, spindex.o2], $
;                   coeffmatrix, ratematrix, heatmatrix

RETURN
END
