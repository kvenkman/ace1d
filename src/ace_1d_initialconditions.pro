; Initial conditions are obtained from an global average atmospheric profile 
; for medium solar activity, at equinox
; The lower boundary O/O2 mixing ratios and TN are initially obtained from MSIS
	restore, 'eqnx_smed.sav'
	restore, 'o_o2_lbc.sav'
	restore, 'tn_lbc.sav'

	nl=0.

	; The values in the save files are on a 29 level pressure grid (TIE-GCM)
	; We interpolate the values to fit them to the ACE1D grid

	; Initial major species (in mmr)
	initial_o = total(reform(o1[nl, 17:18, *, 0]), 1)/2.
	initial_o2 = total(reform(o2[nl, 17:18, *, 0]), 1)/2.
	initial_zp= findgen(29)/2 - 6.75

	; Initial temperatures
	initial_tn = total(reform(tn[nl, 17:18, *, 0]), 1)/2.
	initial_te = total(reform(te[nl, 17:18, *, 0]), 1)/2.
	initial_ti = total(reform(ti[nl, 17:18, *, 0]), 1)/2.

	; Initial O+ and electron densities
	initial_e = total(reform(n_e[nl, 17:18, *, 0]), 1)/2.
	initial_o_p = total(reform(op[nl, 17:18, *, 0]), 1)/2.
	initial_o_p[-1] = initial_e[-1] ; Removing bad values

	; Correcting for top boundary values (as these are initialized to 1e36 in TIEGCM)
	itop=n_elements(initial_tn)

	dtn = initial_tn[-2] - initial_tn[-3]; deriv(initial_zp, initial_tn)
	dte = initial_te[-2] - initial_te[-3]; deriv(initial_zp, initial_te)
	dti = initial_ti[-2] - initial_ti[-3]; deriv(initial_zp, initial_ti)

	initial_tn[itop-1]= initial_tn[-2] + dtn
	initial_te[itop-1]= initial_te[-2] + dte
	initial_ti[itop-1]= initial_ti[-2] + dti

	initial_no = total(reform(no[nl, 17:18, *, 0]), 1)/2.
	initial_n4s = total(reform(n4s[nl, 17:18, *, 0]), 1)/2.
	initial_n2d = total(reform(n2d[nl, 17:18, *, 0]), 1)/2.

	IF(model.nlev NE 29) THEN BEGIN
		zmaj.o = interpol(initial_o, initial_zp, model.zp) > 1e-6
		zmaj.o2 = interpol(initial_o2, initial_zp, model.zp) > 1e-6

		zion.o_p = interpol(initial_o_p, initial_zp, model.zp)>0.0
		zion.e   = interpol(initial_e, initial_zp, model.zp)>0.0
		zminor.no = interpol(initial_no, initial_zp, model.zp)>0.0
		zminor.n4s = interpol(initial_n4s, initial_zp, model.zp)>0.0
		zminor.n2d = interpol(initial_n2d, initial_zp, model.zp)>0.0

		zmaj.tn = interpol(initial_tn, initial_zp, model.zp)>0.0
		zion.te = interpol(initial_te, initial_zp, model.zp)>zmaj.tn
		zion.ti = interpol(initial_ti, initial_zp, model.zp)>zmaj.tn
	ENDIF ELSE BEGIN
		zmaj.o = initial_o
		zmaj.o2 = initial_o2

		zmaj.tn = initial_tn
		zion.o_p = initial_o_p
		zion.e = initial_e
		zion.te = initial_te
		zion.ti = initial_ti
		zminor.no = initial_no
	ENDELSE
		
	; Using MSIS and F10.7 to define the LB mmr for O/O2/N2 and TN
	; The fit coefficients are predetermined from an array of MSIS runs 
	; for different "P" (solar activity index) values
	local_p = 0.5*(inputs.f107d + inputs.f107a)
	o_lb_mmr = o_o2_lbc.o_mmr_fit[0] + local_p*o_o2_lbc.o_mmr_fit[1] + local_p^2*o_o2_lbc.o_mmr_fit[2]
	o2_lb_mmr = o_o2_lbc.o2_mmr_fit[0] + local_p*o_o2_lbc.o2_mmr_fit[1] + local_p^2*o_o2_lbc.o2_mmr_fit[2]
	
	zmaj.o[0] = o_lb_mmr
	zmaj.o2[0] = o2_lb_mmr
	zmaj.n2 = 1. - zmaj.o - zmaj.o2

	zmaj.barm = 1./(zmaj.o/mass.o + zmaj.o2/mass.o2 + zmaj.n2/mass.n2)
	
	zmaj.tn[0] = (tn_lbc.tn_lbc_fit[0] + local_p*tn_lbc.tn_lbc_fit[1] + local_p^2*tn_lbc.tn_lbc_fit[2])
	zion.te = zmaj.tn
	zion.ti = zmaj.tn

	; Converting species from mmr to cm^-3
	zminor.no_mmr = zminor.no
	zminor.no = ace_1d_mmr_cm3(zminor.no , mass.no , zmaj, model, pconst)
	zminor.n4s_mmr = zminor.n4s
	zminor.n4s= ace_1d_mmr_cm3(zminor.n4s, mass.n4s, zmaj, model, pconst)
	zminor.n2d= ace_1d_mmr_cm3(zminor.n2d, mass.n4s, zmaj, model, pconst)

	zmaj.oden  = ace_1d_mmr_cm3(zmaj.o, mass.o, zmaj, model, pconst)
	zmaj.o2den = ace_1d_mmr_cm3(zmaj.o2, mass.o2, zmaj, model, pconst)
	zmaj.n2den = ace_1d_mmr_cm3(1. - zmaj.o - zmaj.o2, mass.n2, zmaj, model, pconst)
	
	; Initialize altitude and gravity values of the pressure grid
    ace_1d_zcalc, model, zmaj, mass, pconst
    
	;; Building a CO2 profile
	co2_lb = 350.e-6 ; Assuming this is the CO2 mmr at the LB of the model
	nco2_lb = co2_lb*zmaj.barm[0]*(zmaj.oden[0] + zmaj.o2den[0] + zmaj.n2den[0])/mass.co2 ; Converting to number density
	zminor.co2[0] 	= nco2_lb
	co2_sch = (pconst.boltz*zmaj.tn)/(mass.co2*pconst.grav/pconst.avo) ; Defining the scale height for CO2
	w = model.dz*(pconst.gask*zmaj.tn)/(pconst.grav*zmaj.barm)
	for i = 1, model.nlev-1 DO zminor.co2[i] = zminor.co2[i-1]*exp(-w[i-1]/co2_sch[i-1])

	; Saving initial conditions in a separate structure for posteriority
	initial = {no:zminor.no, n4s:zminor.n4s, tn:zmaj.tn, te:zion.te, $ 
			   ti:zion.ti, o:zmaj.oden, o2:zmaj.o2den, n2:zmaj.n2den, $
			   e:zion.e, o_p:zion.o_p, zz:zmaj.zz}


