PRO create_initial_conditions
    ; Create files to be used as initial conditions by the ACE1D model
	yyyyddd = 1999080
	latbin = 5.*findgen(37) - 90.
	altbin = findgen(600)+95.
	f107 = 150
	longbin = 15.*findgen(24)

	tn_global_averages = fltarr(n_elements(altbin))
	o1_global_averages = fltarr(n_elements(altbin))
	o2_global_averages = fltarr(n_elements(altbin))
	n2_global_averages = fltarr(n_elements(altbin))
	zp_global_averages = fltarr(n_elements(altbin))

	dlat = 5 * !dtor ; lat bins are 5 degrees apart. Convert to radians
	dlon = 15 * !dtor ; lon bins are 15 degrees apart. Convert to radians
	r_e = 6371. ; Earth radius in km
	coslat = cos(latbin * !dtor) > 0. ; Making sure the cosine goes to zero at +/- 90. 

	ap = 5 ; arbitrarily low geomagnetic activity

	; Integrating and averaging quantities over the earth surface
	; Radius terms gets canceled out
	for j = 0, n_elements(latbin) - 1 do begin
	    for k = 0, n_elements(longbin) - 1 do begin
	        d=msis2k(altbin, yyyyddd, 0., latbin[j], longbin[k], 24.*k/n_elements(longbin), f107, f107, ap)
	        tn_global_averages = tn_global_averages + d.t*coslat[j]*dlat*dlon/(4.*!pi)
	        o1_global_averages = o1_global_averages + d.o*coslat[j]*dlat*dlon/(4.*!pi)
	        o2_global_averages = o2_global_averages + d.o2*coslat[j]*dlat*dlon/(4.*!pi)
	        n2_global_averages = n2_global_averages + d.n2*coslat[j]*dlat*dlon/(4.*!pi)
	    endfor ; loop over longitude
	endfor ; loop over latitude
    
    ; Save outputs in netcdf file
    fid = ncdf_create('initial_conditions.nc', /clobber)
    ncdf_control, fid, /fill
    alt_id = ncdf_dimdef(fid, 'alt', n_elements(altbin))
    
    a_id  = ncdf_vardef(fid, 'alts', [alt_id])
    tn_id = ncdf_vardef(fid, 'tn_global_average', [alt_id])
    o1_id = ncdf_vardef(fid, 'o1_global_average', [alt_id])
    o2_id = ncdf_vardef(fid, 'o2_global_average', [alt_id])
    n2_id = ncdf_vardef(fid, 'n2_global_average', [alt_id])
    
    ncdf_attput, fid, a_id ,  'UNITS', 'km'
    ncdf_attput, fid, tn_id , 'UNITS', 'K'
    ncdf_attput, fid, o1_id , 'UNITS', 'cm^-3'
    ncdf_attput, fid, o2_id , 'UNITS', 'cm^-3'
    ncdf_attput, fid, n2_id , 'UNITS', 'cm^-3'
    
    ncdf_control, fid, /endef
    
    ncdf_varput, fid, a_id, altbin
    ncdf_varput, fid, tn_id,  tn_global_averages
    ncdf_varput, fid, o1_id,  o1_global_averages
    ncdf_varput, fid, o2_id,  o2_global_averages
    ncdf_varput, fid, n2_id,  n2_global_averages
    
    ncdf_close, fid
END
