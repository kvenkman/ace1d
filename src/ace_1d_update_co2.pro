PRO ace_1d_update_co2, zminor, zmaj, model, mass, pconst
; Update CO2 profile based on changes to the background atmosphere

	co2_lb = 350.e-6 ; Assuming this is the CO2 mmr at the bottom of the model
	nco2_lb = co2_lb*zmaj.barm[0]*(zmaj.oden[0] + zmaj.o2den[0] + zmaj.n2den[0])/mass.co2 ; Converting to number density
	zminor.co2[0] 	= nco2_lb
	co2_sch = (pconst.boltz*zmaj.tn)/(mass.co2*pconst.grav/pconst.avo) ; Defining the scale height for CO2
	w = model.dz*(pconst.gask*zmaj.tn)/(pconst.grav*zmaj.barm)
	for i = 1, model.nlev-1 DO zminor.co2[i] = zminor.co2[i-1]*exp(-w[i-1]/co2_sch[i-1])

END
