; Initial conditions for the background neutral atmosphere are obtained from a 
; MSIS based global average atmospheric profile

; For minor species, we assume photochemical equilibrium
; For the ionosphere, we assume photochemical equilibrium, and set the 
; ion and electron temperatures to be equal to the neutral temperature

; The values are interpolated onto a vertical pressure based coordinate system
fid = ncdf_open('../misc/initial_conditions.nc')
tn_id = ncdf_varid(fid, 't_global_average')
o1_id = ncdf_varid(fid, 'o1_global_average')
o2_id = ncdf_varid(fid, 'o2_global_average')
n2_id = ncdf_varid(fid, 'n2_global_average')

ncdf_varget, fid, tn_id, tn
ncdf_varget, fid, o1_id, o1
ncdf_varget, fid, o2_id, o2
ncdf_varget, fid, n2_id, n2

; Total atmospheric number densities
tot_density = o1 + o2 + n2

; Convert number density and temperature to pressure
pressure = tot_density * pconst.boltz * tn

; Interpolating to model coordinates 
p_ref = 5e-4
msis_z = -alog(pressure/p_ref)

zmaj.tn = interpol(tn, msis_z, model.zp)
zmaj.oden = interpol(o1, msis_z, model.zp)
zmaj.o2den = interpol(o2, msis_z, model.zp)
zmaj.n2den = interpol(n2, msis_z, model.zp)

; Establish geometric height of these pressure levels

@ace_1d_minor_species_initial
@ace_1d_ion_species_initial

END
