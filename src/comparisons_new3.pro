;PRO comparisons_new3, output
	inputs = {run_year:1999, start_day:80, ndays:14, runlat:0, timestep:120., save_res:2*24l,f107d:70,f107a:70}
	ace_1d_mainprogram, inputs, output
	
    restore, 'globalave.sav'
    boltz = 1.38e-16
    p_ref = 5e-4
    zp = 0.25*findgen(57) - 7

    ;n_250 = glbmean_250.xno2[*, -1] + glbmean_250.xno[*, -1] + glbmean_250.xnn2[*, -1]
    n_70 = n2_global_averages[* , 0] + o2_global_averages[* , 0] + o_global_averages[* , 0]
    
    ;p_250 = n_250*boltz*glbmean_250.tn[*, -1]
    p_70 = n_70*boltz*t_global_averages[*, 0]

    ;msis_250_zp = -alog(p_250/p_ref)
    msis_70_zp = -alog(p_70/p_ref)

    !p.multi = [0, 2, 2, 0]
    zz = output.zmaj.zz
    @comparisons_new3_plot        
;        print, zmaj.tn
;        print, glbmean_250.tn[*, -1]

        ;print, msis_250_zp
        print, msis_70_zp
END     
