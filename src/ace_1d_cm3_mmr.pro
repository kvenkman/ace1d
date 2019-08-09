FUNCTION ace_1d_cm3_mmr, x, x_mass, zmaj, model, pconst
;  This function converts species from number densities to mass mixing ratios
;  Required parameters are
; - number density of species to be converted, along with molecular mass
; - mean molecular mass
; - neutral temperature
;
; The function will return an array containing mass mixing ratios
;  with dimensions same as the input array

; The returned array will be dimensionless

;"mmr" is mass mixing ratio with (o2+o+n2=1). The formula from mmr to convert to cm-3 is: f(cm3) = f(mmr) * pkt * barm / w
; where pkt = p0*e(-z)/kT (k is boltzman's, T is neutral temperature on the history), barm is the mean mass,
; and w is the molecular weight of the species being converted. barm is calculated as 1./(o2/32+o/16+n2/28), 
; where o2,o,n2 are in mmr.

pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmaj.tn)
x_mmr = x*x_mass/(zmaj.barm*pkt)

return, x_mmr ; 
END
