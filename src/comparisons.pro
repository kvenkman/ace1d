restore, 'globalave.sav'
;restore, 'model_output_tn_85gw.sav'
restore, 'model_output_tn_70gw.sav'
;restore, 'model_output_tn_70gw_eddy_test.sav'
;restore, 'model_output_tn_70gw_const_tn.sav'
;restore,'model_output_tn_70gw_no_minorsolv.sav'
;restore, 'model_output_tn_70gw_halfo3p.sav'
;restore, 'model_output_tn_70gw_zero_nov.sav'
;restore, 'model_output_tn_70gw_half_o3pcooling.sav'
;restore, 'model_output_tn_85gw.sav'
;restore, 'model_output_tn_90gw.sav'
;restore, 'model_output_tn_70gw_withoutnov.sav'
;restore, 'glbmean_tntest.sav' ; Contains GMM TN
restore, 'tiegcm_f107_tn.sav'

first = 0
last = -1
!p.multi = [0, 2, 2, 0]
n_window = 0
scale = 1
window, n_window+1, xsize = scale*1024, ysize = scale*768
cgplot, alog10(n2_densities[* ,first]), zz[* ,0],xr=[3, 14], title = 'Major species, F10.7=70',xtitle='Log Densities (cm!E-3!N)',ytitle = 'Altitude', yr = [100, 400]
cgoplot, alog10(o2_densities[* ,first]), zz[* ,0], color = 'blue'
cgoplot, alog10(o_densities[* ,first]), zz[* ,0], color = 'red'

cgoplot, alog10(n2_global_averages[* ,first]), altbin, linestyle = 2
cgoplot, alog10(o2_global_averages[* ,first]), altbin, linestyle = 2, color = 'blue'
cgoplot, alog10(o_global_averages[* ,first]), altbin, linestyle = 2, color = 'red'
plots, [10, 11.5], [350, 350]
plots, [10, 11.5], [300, 300], linestyle = 2
xyouts, 12, 350, 'ACE'
xyouts, 12, 300, 'MSIS'


cgplot, alog10(n2_densities[* ,last]), zz[* ,last], title = 'Major species, F10.7=250',xtitle='Log Densities (cm!E-3!N)',ytitle = 'Altitude', yr = [100, 550]
cgoplot, alog10(o2_densities[* ,last]), zz[* ,last], color = 'blue'
cgoplot, alog10(o_densities[* ,last]), zz[* ,last], color = 'red'
xyouts, 12, 500, 'O',color=11
xyouts, 12, 450, 'O2',color = 3
xyouts, 12, 400, 'N2'

cgoplot, alog10(n2_global_averages[* ,last]), altbin, linestyle = 2
cgoplot, alog10(o2_global_averages[* ,last]), altbin, linestyle = 2, color = 'blue'
cgoplot, alog10(o_global_averages[* ,last]), altbin, linestyle = 2, color = 'red'

cgplot, t_global_averages[*, first], altbin, linestyle = 2, yr = [100, 200], ytitle = 'Altitude (kms)', xtitle = '(K)',title='TN, F10.7=70'
cgoplot, tn[*, first], zz[*, first]

cgplot, t_global_averages[*, last], altbin, linestyle = 2, yr = [100, 200], xr = [0, 1500], ytitle = 'Altitude (kms)', xtitle = '(K)',title='TN, F10.7=250'
cgoplot, tn[*, last], zz[*, last]
;tvpng,'plot1.png'
print, "MSIS Exo : ", t_inf_averages
print, "Model Exo : ", t_exo

!p.multi = [0, 2, 1, 0]
window, n_window+2, xsize = scale*1024, ysize = scale*384
cgplot, no_densities[*, first]/1e7, zz[*, 0], yr = [100, 200], /xs, xr = [0, 25], title = 'NO Densities', ytitle = 'Altitude (kms)',color = 'blue',xtitle = '(x1E7 cm!E-3!N)'
cgoplot, no_densities[*, last]/1e7, zz[*, last], yr = [100, 200], /xs,color = 'red'
plots, [0,25], [110, 110], linestyle = 2
xyouts,0.5, 111, '110 km'
xyouts,8, 180, 'F10.7 = 70', color = 3
xyouts,8, 165, 'F10.7 = 250', color = 11

cgplot, f107, tn[-1, *], xr  = [60, 260], color = 'blue',xtitle='F10.7',ytitle='T!DEXO!N (K)',title='T!DEXO!N vs F10.7', yr = [600, 1600]
cgoplot, f107, t_global_averages[-1, *], psym = 4
;cgoplot, f107, glb_tn, color = 'red'
cgoplot, f107, tiegcm_tn[-1, *], color = 'red'

xyouts, 100, 1500, 'TIEGCM', color=11
;xyouts, 200, 300, 'GMM', color=11
xyouts, 100, 1300, 'ACE',color = 3
xyouts, 100, 1400, 'MSIS'
;tvpng,'plot2.png'

!p.multi = 0
