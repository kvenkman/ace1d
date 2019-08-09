    cgplot, alog10(zmaj.n2den), model.zp, yr = [-7,5], xr = [3,14]
    cgoplot, alog10(zmaj.o2den), model.zp, color = 'blue'
    cgoplot, alog10(zmaj.oden), model.zp, color = 'red'
    cgoplot, alog10(initial.n2), model.zp, linestyle=2 
    cgoplot, alog10(initial.o2), model.zp, color = 'blue', linestyle=2 
    cgoplot, alog10(initial.o), model.zp, color = 'red', linestyle=2 
    
    cgplot, zion.te, model.zp, color = 'red', xr = [100, 2500], yr = [-7, 5]
    cgoplot, zion.ti, model.zp, color = 'blue'
    cgoplot, zmaj.tn, model.zp
    cgoplot, initial.tn, model.zp, linestyle = 2
    cgoplot, initial.te, model.zp, color = 'red', linestyle = 2
    cgoplot, initial.ti, model.zp, color = 'blue', linestyle = 2
    
    cgplot, alog10(zion.e), model.zp, xr = [2,7],yr = [-7,5], linestyle = 2,title = 'Ion species',xtitle='Log!D10!N densities',ytitle = 'Pressure level'
    cgoplot, alog10(zion.o_p), model.zp, color = 'purple'
    cgoplot, alog10(zion.n_p), model.zp, color = 'blue'
    cgoplot, alog10(zion.no_p), model.zp, color = 'red'
    cgoplot, alog10(zion.n2_p), model.zp, color = 'green'
    cgoplot, alog10(zion.o2_p), model.zp, color = 'orange'

    cgplot, alog10(zminor.no), model.zp, xr = [0, 8], yr = [-7, 5],title = 'Minor species',xtitle='Log!D10!N densities',ytitle = 'Pressure level'
    cgoplot, alog10(zminor.n2d), model.zp, color = 'red'
    cgoplot, alog10(zminor.n4s), model.zp, color = 'blue'
    cgoplot, alog10(zminor.n2p), model.zp, color = 'orange'
    cgoplot, alog10(zminor.o1d), model.zp, color = 'red', psym = 4
