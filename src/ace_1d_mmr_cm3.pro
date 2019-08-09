FUNCTION ace_1d_mmr_cm3, x, x_mass, zmaj, model, pconst

; This function converts species from mass mixing ratios (mmr) to number densities
;  Required parameters are
; - mmr of species to be converted, along with molecular mass
; - mmr of O and O2
; - neutral temperature
;
; The function will return an array containing number densities with dimensions same as the input array
; Units on the returned array will be #/cm3

;"mmr" is mass mixing ratio with (o2+o+n2=1). The formula from mmr to convert to cm-3 is: f(cm3) = f(mmr) * pkt * barm / w
;where pkt = p0*e(-z)/kT (k is boltzman's, T is neutral temperature on the history), barm is the mean mass,
;and w is the molecular weight of the species being converted. barm is calculated as 1./(o2/32+o/16+n2/28), 
;where o2,o,n2 are in mmr.

; The input is a single column


x_cm3 = x*zmaj.barm*model.p0*exp(-model.zp)/(zmaj.tn*pconst.boltz*x_mass)
return, x_cm3 ; 
END
