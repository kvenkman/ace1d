PRO ace_1d_zcalc, model, zmaj, mass, pconst
	; Calculate altitudes of pressure levels

	zmaj.z[0] = model.zbound
	pconst.grav[0] = pconst.gravref*(pconst.re/(pconst.re+zmaj.z[0]))^2.

	; Adding column height at each vertical grid point to obtain the altitude of the next grid point
	FOR i = 1, model.nlev - 1 DO BEGIN
	  ;w = model.dz*(pconst.gask*zmaj.tn[i-1])/(pconst.grav*zmaj.barm[i-1])
	  w = model.dz*(pconst.gask*zmaj.tn[i-1])/(pconst.grav[i-1]*zmaj.barm[i-1])
	  zmaj.z[i] = zmaj.z[i-1] + w
	  pconst.grav[i] = pconst.gravref*(pconst.re/(pconst.re+zmaj.z[i]))^2.
	ENDFOR

	zmaj.zz=zmaj.z/1e5 ; km
END
