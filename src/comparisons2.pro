restore, 'tiegcm_global_avg.sav'
restore, 'globalave.sav'
;	restore, 'glbmean_majors_80.sav'
restore, 'glbmean_majors_250.sav'

; contents of glbmean save
; xno2, xnn2, xno, xno1d, zpht, xom, xo2m, xn2m, xnoc, qjoule, tn, te, ti, tn_msis, lev,

!p.multi = [0, 2, 2, 0]
last = -1

cgplot, alog10(zmaj.n2den), zmaj.zz, title = 'Major species',xtitle='Densities (cm!E-3!N)',ytitle = 'Altitude (kms)'
cgoplot, alog10(zmaj.o2den), zmaj.zz, color = 'blue'
cgoplot, alog10(zmaj.oden), zmaj.zz, color = 'red'

;cgoplot, alog10(o1_ave), zz_ave, color = 'red', linestyle = 2
;cgoplot, alog10(o2_ave), zz_ave, color = 'blue', linestyle = 2
;cgoplot, alog10(n2_ave), zz_ave, linestyle = 2

cgoplot, alog10(xno[*, last]), zpht[*, last], color = 'red', linestyle = 2
;cgoplot, alog10(xno2[*, last]), zpht[*, last], color = 'blue', linestyle = 2
;cgoplot, alog10(xnn2[*, last]), zpht[*, last], linestyle = 2
cgoplot, alog10(xo2m), zpmht, color = 'blue', linestyle = 2
cgoplot, alog10(xn2m), zpmht, linestyle = 2 ; MSIS from glbmean. Just to check.

cgoplot, alog10(n2_global_averages[* ,last]), altbin, psym = 4
cgoplot, alog10(o2_global_averages[* ,last]), altbin, psym = 4, color = 'blue'
cgoplot, alog10(o_global_averages[* ,last]), altbin, psym = 4, color = 'red'


cgplot, zmaj.tn, zmaj.zz, xr = [0, 1450], title = 'Neutral Temperature (K)', yr = [100, 500]
;cgoplot, tn_ave, zz_ave, linestyle = 2
;cgoplot, tn[*, last], zpht[*, last], linestyle = 2
cgoplot, tn_msis, zpmht, linestyle = 2 ; MSIS from glbmean. Just to check.
cgoplot, t_global_averages[*, last], altbin, psym = 4

cgplot, alog10(zmaj.n2den), model.zp, title = 'Major species',xtitle='Densities (cm!E-3!N)',ytitle = 'Pressure level'
cgoplot, alog10(zmaj.o2den), model.zp, color = 'blue'
cgoplot, alog10(zmaj.oden), model.zp, color = 'red'

;cgoplot, alog10(o1_ave), lev, color = 'red', linestyle = 2
;cgoplot, alog10(o2_ave), lev, color = 'blue', linestyle = 2
;cgoplot, alog10(n2_ave), lev, linestyle = 2

cgoplot, alog10(xno[*, last]), lev, color = 'red', linestyle = 2
cgoplot, alog10(xno2[*, last]), lev, color = 'blue', linestyle = 2
cgoplot, alog10(xnn2[*, last]), lev, linestyle = 2

cgplot, zmaj.tn, model.zp, xr = [0, 1450], title = 'Neutral Temperature (K)'
cgoplot, tn[*, last], lev, linestyle = 2
;cgoplot, tn_ave, lev, linestyle = 2

!p.multi = 0
