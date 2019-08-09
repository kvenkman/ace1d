pro modatm,atm,season=season

; procedure to read model atmosphere files
; us standard atmosphere is nominal case
; if season is supplied, then options are:
; 'tro' for tropical
; 'mls' for midlatitude summer
; 'mlw' for midlatitude winter
; 'sas' for subarctic summer
; 'saw' for subarctic winter
; see www.atm.ox.ac.uk/RFM/atm

; results are returned as a structure

file='std.atm'
if keyword_set(season) then file=strtrim(season,2)+'.atm'

n_text=3
n_alt=50
s='a string'

z=fltarr(n_alt) & t=fltarr(n_alt) & p=fltarr(n_alt) & h2o=fltarr(n_alt) & co2=fltarr(n_alt) & o3=fltarr(n_alt)
n2o=fltarr(n_alt) & co=fltarr(n_alt) & ch4=fltarr(n_alt) & o2=fltarr(n_alt)

openr,lun,file,/get_lun
for i=1,n_text do readf,lun,s

readf,lun,s
readf,lun,z
readf,lun,s
readf,lun,p
readf,lun,s
readf,lun,t
readf,lun,s
readf,lun,h2o
readf,lun,s
readf,lun,co2
readf,lun,s
readf,lun,o3
readf,lun,s
readf,lun,n2o
readf,lun,s
readf,lun,co
readf,lun,s
readf,lun,ch4
readf,lun,s
readf,lun,o2

ar=1e4+fltarr(n_alt)

tot=(o2+ch4+co+n2o+o3+co2+h2o+ar) / 1e6 ; converts from parts per million to fraction
n2=(1.0 - tot)*1e6 ; N2 density is the rest of the atmosphere

mxr = {n2:n2,o2:o2,ar:ar,ch4:ch4,co:co,n2o:n2o,o3:o3,co2:co2,h2o:h2o}

; now we need to get number densities
; use the ideal gas law to convert p&t to n
; p is in mb so conver to pascals 
n=(p/1000.*!a.conv.pascal_per_bar) / (!a.const.kb*t) ; n in m^-3
n=n/1e6 ; n in cm^-3


scale = n/1e6 ; this is scale factor so mixing ratio in parts per million is converted to number density
conc = {n2:n2*scale,o2:o2*scale,ar:ar*scale,ch4:ch4*scale,co:co*scale,n2o:n2o*scale,o3:o3*scale,co2:co2*scale,h2o:h2o*scale}

; calculate column densities

dz=z(1:*)-z
dz=[dz,dz(n_alt-2)]*1e5

cden_n=fltarr(n_alt) & cden_n2=fltarr(n_alt) & cden_o2=fltarr(n_alt) & cden_ar=fltarr(n_alt) & cden_co=fltarr(n_alt)
cden_n2o=fltarr(n_alt) & cden_o3=fltarr(n_alt) & cden_co2=fltarr(n_alt) & cden_h2o=fltarr(n_alt) & cden_ch4=fltarr(n_alt)

scale = n[n_alt-1]/1e6
cden_n[n_alt-1] = n[n_alt-1]*dz[n_alt-1]*scale
cden_n2[n_alt-1] = n2[n_alt-1]*dz[n_alt-1]*scale
cden_o2[n_alt-1] = o2[n_alt-1]*dz[n_alt-1]*scale
cden_ar[n_alt-1] = ar[n_alt-1]*dz[n_alt-1]*scale
cden_co[n_alt-1] = co[n_alt-1]*dz[n_alt-1]*scale
cden_n2o[n_alt-1] = n2o[n_alt-1]*dz[n_alt-1]*scale
cden_o3[n_alt-1] = o3[n_alt-1]*dz[n_alt-1]*scale
cden_co2[n_alt-1] = co2[n_alt-1]*dz[n_alt-1]*scale
cden_h2o[n_alt-1] = h2o[n_alt-1]*dz[n_alt-1]*scale

for i=n_alt-2,0,-1 do begin
scale = n[i]/1e6
cden_n[i] = cden_n[i+1] + n[i]*dz[i]*scale
cden_n2[i] = cden_n2[i+1] + n2[i]*dz[i]*scale
cden_o2[i] = cden_o2[i+1] + o2[i]*dz[i]*scale
cden_ar[i] = cden_ar[i+1] + ar[i]*dz[i]*scale
cden_ch4[i] = cden_ch4[i+1] + ch4[i]*dz[i]*scale
cden_co[i] = cden_co[i+1] + co[i]*dz[i]*scale
cden_n2o[i] = cden_n2o[i+1] + n2o[i]*dz[i]*scale
cden_o3[i] = cden_o3[i+1] + o3[i]*dz[i]*scale
cden_co2[i] = cden_co2[i+1] + co2[i]*dz[i]*scale
cden_h2o[i] = cden_h2o[i+1] + h2o[i]*dz[i]*scale
endfor

vcden={n:cden_n,n2:cden_n2,o2:cden_o2,ar:cden_ar,ch4:cden_ch4,co:cden_co,n2o:cden_n2o,o3:cden_o3,co2:cden_co2,h2o:cden_co2}
atm = {n_alt:n_alt,z:z,dz:dz,n:n,t:t,p:p,conc:conc,mxr:mxr,vcden:vcden}

return
end





