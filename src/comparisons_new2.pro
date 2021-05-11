;PRO comparisons_new2, output
	inputs = {run_year:1999, start_day:80, ndays:14, runlat:0, timestep:120., save_res:2*24l,f107d:250,f107a:250}
	ace_1d_mainprogram, inputs, output
	
    restore, 'globalave_test.sav'
    boltz = 1.38e-16
    p_ref = 5e-4
    zp = 0.25*findgen(57) - 7

    ;n_250 = glbmean_250.xno2[*, -1] + glbmean_250.xno[*, -1] + glbmean_250.xnn2[*, -1]
    n_250 = n2_global_averages[* ,-1] + o2_global_averages[* ,-1] + o_global_averages[* ,-1]
    
    ;p_250 = n_250*boltz*glbmean_250.tn[*, -1]
    p_250 = n_250*boltz*t_global_averages[*, -1]

    msis_250_zp = -alog(p_250/p_ref)
    
     ;cgplot, output.zmaj.oden, output.zmaj.zz, color = 'blue', yr = [90, 200]
     cgplot, output.azmajor[10].oden, zp, yr = [-8, 8], /xl;, xr = [1e1, 1e14]
     cgoplot, output.azmajor[-10].oden, zp, color = 'blue'

     ;cgoplot, output.azmajor[10].o2den, zp, linestyle = 2
     ;cgoplot, output.azmajor[-10].o2den, zp, color = 'blue', linestyle = 2
    
    ;STOP
    
    zz = output.zmaj.zz

    @comparisons_new2_plot
;        print, zmaj.tn
;        print, glbmean_250.tn[*, -1]
        print, msis_250_zp
        STOP
END     
