; Initial conditions for the background neutral atmosphere are obtained from a 
; MSIS based global average atmospheric profile

; For minor species, we assume photochemical equilibrium
; For the ionosphere, we assume photochemical equilibrium, and set the 
; ion and electron temperatures to be equal to the neutral temperature

; The values are interpolated onto a vertical pressure based coordinate system
; Model constants
model  = {  $
            nlev     : 57.                  ,$ ; 57 vertical grid points
            ztop     : 7.                   ,$ ; +7 is the top-most vertical grid point (pressure level)
            zbot     :-7.                   ,$ ; -7 is the bottom-most vertical grid point
            p0       : 5.e-4                ,$ ; Model reference pressure in dyne/cm^2, corresponding to the pressure at the grid point z = 0          
            zp       : findgen(57)/4 - 7.   ,$ ; The vertical grid has 4 grid points per scale height. 
            zbound   : 97e5                 ,$ ; The lower pressure boundary of the model is assumed to be at 97 km
            dz       : 0.25                 ,$ ; dz = 0.25

	         nday_hours      : 24, $
	         nhour_minutes   : 60, $
	         nminute_seconds : 60  $
         }  
         
; Physical constants and conversion factors
pconst = {  boltz : 1.38066e-16,$       ; boltzmann constant  in ergs K^-1        
            h     : 6.626e-27  ,$       ; in ergs s
            c     : 2.9979e10  ,$       ; in cm s^-1
            avo   : 6.023e23   ,$       ; avogadro's number           
            gask  : 8.314e7    ,$       ; gas constant                
            grav  : fltarr(model.nlev),$  ; acceleration due to gravity varies as a function of height, in cm s^-2
            gravref: 980.      ,$       ; Reference value of gravity
            re    : 6371.e5    ,$       ; Earth radius in cm
            cgm   : 0.5        ,$       
            eff   : 0.33       ,$       ; Solar EUV heating efficiency, from TIE-GCM
            ev2erg: 1.602E-12  ,$		; eV to ergs
            amu   : 1.6605e-24 ,$
            q_e   : 1.602e-19  ,$       ; Electron charge
            m_e   : 9.109e-28   $		; Electron mass
         }      

zmaj = {    o      : fltarr(model.nlev)   ,$
            o2     : fltarr(model.nlev)   ,$
            n2     : fltarr(model.nlev)   ,$
            oden   : fltarr(model.nlev)   ,$
            o2den  : fltarr(model.nlev)   ,$
            n2den  : fltarr(model.nlev)   ,$
            tn     : fltarr(model.nlev)   ,$
            z      : fltarr(model.nlev)   ,$ ; Altitude of pressure levels
            zz     : fltarr(model.nlev)   ,$ ; Altitude in km
            barm   : fltarr(model.nlev)   ,$ ; mean mass
            cp     : fltarr(model.nlev)   ,$
            kt     : fltarr(model.nlev)    $
       }            
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

@ace_1d_minor_species_initial
@ace_1d_ion_species_initial

END
