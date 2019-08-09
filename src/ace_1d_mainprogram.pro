PRO ace_1d_mainprogram, inputs, output
    ; Authors
    ; K Venkataramani
    ; S M Bailey
    ; Center for Space Science and Engineering Research, Virginia Tech
    ; Blacksburg, VA

    ; The present model builds upon the 1D Nitric Oxide model
    ; and extends it to include a solver for
    ; Major Species densities
    ; Other minor species
    ; Ion densities, and
    ; Neutral, ion and electron temperatures


	; Number of timesteps required for the run
	nday_hours = 24
	nhour_minutes = 60
	nminute_seconds = 60

    nsteps = long(inputs.ndays*nday_hours*nhour_minutes/ $
    	         (inputs.timestep/nminute_seconds) + 1)

	; Number of saves in the output
    n_saves = ((long(nsteps)) / inputs.save_res)+1

	; Initialize constants & variable structures
    @ace_1d_model_parameters
    @ace_1d_loadeuv
    @ace_1d_defvariables
    @ace_1d_chem_parameters
    @ace_1d_initialconditions

	;;;;;;
	; These parameters do not currently affect model run
    model_time.day = inputs.start_day
    model_time.year = inputs.run_year
    model_time.time = 43200.00 ; Run starts at noon
    ;;;;;;

	; Set solar activity levels
	model_sun.f107d = inputs.f107d
    model_sun.f107a = inputs.f107a

    ; Initialize solar spectrum
    ace_1d_euvac, solspec, model_sun.f107d, model_sun.f107a, zpid, zcol, inputs

    ; Initialize model chemistry
    ace_1d_updatechem, zmaj, zion, zminor, model, $
                       spindex, exomatrix, coeffmatrix, ratematrix, $
                       heatmatrix

    initial.zz=zmaj.zz

	; Initialize molecular diffusion and thermal conductivity coefficients 
	; and heat capacity
    ace_1d_km, zmaj, zminor, mass, model, pconst, k_m, flux
    ace_1d_cpkt, zmaj, mass, pconst, model

    ; Counter to monitor history saves
    n_save=0l

;	@ace_1d_initialconditions_loop
    @ace_1d_initialconditions


	FOR i = 1., nsteps - 1 DO BEGIN
	
		ace_1d_chapman, model_sun, zcol, pconst, zmaj, model, mass, i

		ace_1d_solarphotonproc, zmaj, model, pconst, solspec, zcol, zpid, zei, $
		                        zion, heattermse, edep, branching, xsec, $ 
		                        model_sun, i

		ace_1d_nophot,model,model_sun,zmaj,zcol,zpid

		ace_1d_updatechem, zmaj, zion, zminor, model, $
		                   spindex, exomatrix, coeffmatrix, ratematrix, $
		                   heatmatrix

		; the "now" variables are the values for this time step
		; the values for the next time step are calculated 
		; in the following routines
		zmajnow   = zmaj
		zminornow = zminor
		zionnow   = zion

		ace_1d_chemistry_ions, zmajnow, zminornow, zionnow, zion,$
		                       ratematrix, coeffmatrix, yieldmatrix, $
		                       spindex, zpid,zei,model,pconst,mass, i

		ace_1d_chemistry_oddn, zmajnow, zminornow, zionnow, zminor, zpid, zei, $
		                       ratematrix, coeffmatrix, yieldmatrix, model, $ 
		                       spindex, termsn2p, termsn2d, termsn4s, termsno, $
		                       flux, k_m, mass, pconst, model_sun, i

		ace_1d_chemistry_ox, zmajnow, zminornow, zionnow, zminor, zpid,  $
							 coeffmatrix, ratematrix, yieldmatrix, branching, $ 
							 spindex, model

		ace_1d_n2a, zminornow, zmajnow, zei, zminor, coeffmatrix, spindex, model

		ace_1d_update_co2, zminor, zmaj, model, mass, pconst

		; chem heating, euvheating, eionheating must be in that order....
		ace_1d_chem_heating, heatmatrix, spindex, heatterms, model, i
		
		ace_1d_euvheating, zmajnow, zcol, xsec, branching, solspec, heatterms, $
						   pconst, model, model_sun
		
		ace_1d_jouleheating, zmajnow, zminornow, zionnow, pconst, mass, model, $
							 heatterms, model_sun, i
		
		ace_1d_auroralheating, zmajnow, heatterms, model, pconst
		ace_1d_cooling, zmajnow, zminornow, zionnow, termsno, pconst, model, $
						coolterms

		ace_1d_eionheating, zmajnow, zminornow, zionnow, $
		                    heattermse, cooltermse, heattermsi, cooltermsi, $ 
		                    heatterms, heatmatrix, spindex, pconst, model, i

		ace_1d_major_chem, zmajnow, zminornow, zionnow, zei, zpid, coeffmatrix,$
		                   yieldmatrix, spindex, o_o2_terms, i, mass, model, $
		                   pconst

		ace_1d_te_ti, model, pconst, zion, zmajnow, zionnow, $ 
					  heattermse.q_total, cooltermse, heattermsi, cooltermsi, $ 
					  heatterms
					  
		ace_1d_dt, zmaj, zmajnow, zionnow, model, pconst, mass, heatterms, $ 
				   coolterms, cooltermse, cooltermsi

		ace_1d_major_solv, zmajnow, zmaj, pconst, model, mass, o_o2_terms, i

		ace_1d_zcalc, model, zmaj, mass, pconst

		ace_1d_km, zmaj, zminor, mass, model, pconst, k_m, flux

		ace_1d_cpkt, zmaj, mass, pconst, model

		; Set up the output variables on the first iteration of the model run
		IF i eq 1 THEN BEGIN
		    aflux         =replicate(flux      ,n_saves)
		    acoolterms    =replicate(coolterms ,n_saves)
		    aheatterms    =replicate(heatterms ,n_saves)
		    acooltermse    =replicate(cooltermse ,n_saves)
		    azion    =replicate(zion           ,n_saves)
		    azminor  =replicate(zminor         ,n_saves)
		    azmajor  =replicate(zmaj           ,n_saves)
		    azpid    = replicate(zpid           ,n_saves)
		    azei     = replicate(zei            ,n_saves)
		    azcol    = replicate(zcol           ,n_saves)
		    ao_o2_terms = replicate(o_o2_terms  ,n_saves)
		ENDIF

		; Save model variables at the specified time resolution
		IF long(i-1) mod long(inputs.save_res) eq 0 THEN BEGIN
		    aflux[n_save] = flux
		    aheatterms[n_save] = heatterms
		    acoolterms[n_save] = coolterms
		    acooltermse[n_save] = cooltermse
		    azion[n_save]    =zion
		    azminor[n_save]  =zminor
		    azmajor[n_save]  =zmaj
		    azpid[n_save] = zpid
		    azei[n_save] = zei
		    azcol[n_save] = zcol

			n_save=n_save+1l
		ENDIF

	; END OF MAIN LOOP
	ENDFOR

    ; This is the final output of the model
	output={zmaj        : zmaj       , $
			zminor      : zminor     , $
			zion        : zion       , $
			coeffmatrix : coeffmatrix, $
			coolterms   : coolterms  , $
			cooltermse  : cooltermse , $
			cooltermsi  : cooltermsi , $
		    heatterms   : heatterms  , $
		    heattermse  : heattermse , $
		    heatmatrix  : heatmatrix , $
		    zpid        : zpid       , $
		    zei         : zei        , $
		    zcol        : zcol       , $
		    spindex     : spindex    , $
		    yieldmatrix : yieldmatrix  $
		    }
END
