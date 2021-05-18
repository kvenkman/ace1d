	!p.multi = [0, 2, 2, 0]

	test_altbin = 90. + 10.*findgen(70)
	zp_test = interpol(zp, zz, test_altbin)
	
	myytickv = [zp_test[1], zp_test[3], zp_test[6], zp_test[11], zp_test[21], zp_test[31], zp_test[41], zp_test[51]]
	myytickname = strtrim(floor([test_altbin[1],test_altbin[3],test_altbin[6],test_altbin[11],test_altbin[21],test_altbin[31], $
				test_altbin[41],test_altbin[51]]), 2)

    window, 1, xsize = 1024, ysize = 2.*384
        cgplot, output.zmaj.n2den, zp, /xl, yr = [-7, 5], title = 'F10.7 = 70', ys = 9
        ;cgoplot, glbmean_250.xnn2[*, -1], msis_250_zp, linestyle = 2
        cgoplot, n2_global_averages[*, 0], msis_70_zp, linestyle = 2
        
        cgoplot, output.zmaj.o2den, zp, /xl, yr = [-7, 5], color = 'blue'
        ;cgoplot, glbmean_250.xno2[*, -1], msis_250_zp, linestyle = 2, color = 'blue'
        cgoplot, o2_global_averages[*, 0], msis_70_zp, linestyle = 2, color = 'blue'
        
        cgoplot, output.zmaj.oden, zp, /xl, yr = [-7, 5], color = 'red'
        ;cgoplot, glbmean_250.xno[*, -1], msis_250_zp, linestyle = 2, color = 'red'
        cgoplot, o_global_averages[*, 0], msis_70_zp, linestyle = 2, color = 'red'
		cgaxis, yaxis = 1, ytickv = myytickv, ytickname = myytickname, /ys, yminor = 1, yticks = 7        

        cgplot, output.zmaj.n2den, output.zmaj.zz, /xl, yr = [100, 400], ys = 9
        cgoplot, n2_global_averages[*, 0], altbin, linestyle = 2
        
        cgoplot, output.zmaj.o2den, output.zmaj.zz, /xl, yr = [-7, 5], color = 'blue'
        cgoplot, o2_global_averages[*, 0], altbin, linestyle = 2, color = 'blue'

        cgoplot, output.zmaj.oden, output.zmaj.zz, /xl, yr = [-7, 5], color = 'red'
        cgoplot, o_global_averages[*, 0], altbin, linestyle = 2, color = 'red'
		cgaxis, yaxis = 1, ytickv = myytickv, ytickname = myytickname, /ys, yminor = 1, yticks = 7

;       cgplot, output.zmaj.tn, zp, yr = [-7, 5], xr = [100, 800]
;       cgoplot, t_global_averages[*, 0], msis_70_zp, linestyle = 2

        ;cgplot, output.zmaj.tn, output.zmaj.zz, yr = [100, 500], xr = [100, 800]
        ;cgoplot, t_global_averages[*, 0], altbin, linestyle = 2
        
        
;       cgplot, output.zmaj.tn, zp, yr = [-7, 7], xr = [100, 2100], ys = 9
;       cgoplot, output.zion.te, zp, color = 'red'
;       cgoplot, output.zion.ti, zp, color = 'blue'
;       cgaxis, /yaxis, /yl, yr = [zz[0], zz[-1]], ytickv = [100, 120, 150, 200, 250, 300, 350, 400, 500], /ys, yticks = 9

		zz = output.zmaj.zz
		
        cgplot, output.zmaj.tn, zp, yr = [-7, 7], xr = [100, 2100], ys = 9
        cgoplot, output.zion.te, zp, color = 'red'
        cgoplot, output.zion.ti, zp, color = 'blue'
        cgaxis, yaxis = 1, ytickv = myytickv, ytickname = myytickname, /ys, yminor = 1, yticks = 7        
        
        cgplot, output.zmaj.tn, zp, yr = [-7, 7], ys = 9, xr = [0, 1000]
        cgoplot, t_global_averages[*, 0], msis_250_zp, linestyle = 2
        cgaxis, yaxis = 1, ytickv = myytickv, ytickname = myytickname, /ys, yminor = 1, yticks = 7        
	 !p.multi = 0