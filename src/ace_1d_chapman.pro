PRO ace_1d_chapman, model_sun, zcol, pconst, zmaj, model, mass, count
; calculate vertical & slant column densities for current zenith angle
; Used for optical depth calculations in ace_1d_solarphotoproc



model_sun.sza =  60.
col_o = fltarr(model.nlev)
col_o2 = fltarr(model.nlev)
col_n2 = fltarr(model.nlev)

col_s_o = fltarr(model.nlev)
col_s_o2 = fltarr(model.nlev)
col_s_n2 = fltarr(model.nlev)

zcol.w = model.dz*(pconst.gask*zmaj.tn)/(pconst.grav*zmaj.barm) ; (kT/mg)*dz (dz = 0.25 scale height)
w = zcol.w

; vertical column densities 
o=zmaj.oden
o2=zmaj.o2den
n2=zmaj.n2den

col_o[-1]  = zmaj.oden[-1]*w[-1]
col_o2[-1] = zmaj.o2den[-1]*w[-1]
col_n2[-1] = zmaj.n2den[-1]*w[-1]
FOR i = model.nlev-2, 0., -1. DO BEGIN
    col_o[i]  = col_o[i+1]  +  zmaj.oden[i]*w[i]
    col_o2[i] = col_o2[i+1] + zmaj.o2den[i]*w[i]
    col_n2[i] = col_n2[i+1] + zmaj.n2den[i]*w[i]
ENDFOR
zcol.o  = col_o
zcol.o2 = col_o2
zcol.n2 = col_n2

; Refer to Solomon & Brasseur, Ch. 4, Box 4.1
term1 = pconst.re + zmaj.z[1:-1]
term2 = pconst.re + zmaj.z[0:-2]

;sin_chi = sin(75. * !dtor)
sin_chi = sin(model_sun.sza * !dtor)
mu = sqrt(term1^2 - (term2^2 * sin_chi^2)) / term1

mu_top = sqrt((pconst.re + zmaj.z[-1] + w[-1])^2. - (sin_chi^2)*(pconst.re + zmaj.z[-1])^2.)/$
         (pconst.re + zmaj.z[-1] + w[-1])

mu = [mu, mu_top]
col_s_o[-1]  =  o[-1]*w[-1]/mu[-1]
col_s_o2[-1] = o2[-1]*w[-1]/mu[-1]
col_s_n2[-1] = n2[-1]*w[-1]/mu[-1]
FOR i = model.nlev-2, 0., -1. DO BEGIN
    col_s_o[i]  = col_s_o[i+1]  +  o[i]*w[i]/mu[i]
    col_s_o2[i] = col_s_o2[i+1] + o2[i]*w[i]/mu[i]
    col_s_n2[i] = col_s_n2[i+1] + n2[i]*w[i]/mu[i]
ENDFOR
zcol.so  = col_s_o
zcol.so2 = col_s_o2
zcol.sn2 = col_s_n2

;zcol.o = zcol.so
;zcol.o2 = zcol.so2
;zcol.n2 = zcol.sn2
END
