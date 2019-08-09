cgplot, zminor.n_4s_pce,zmaj.zz, /xl, /ys,xr = [1e4,1e10], color = 10, yr = [80, 200]
cgoplot, zminor.no_pce,zmaj.zz, color=2
cgoplot, initial.no, zmaj.zz,color = 2, linestyle = 1

cgoplot, zminor.no,zmaj.zz, color=2, linestyle = 2
cgoplot, zminor.n_4s,zmaj.zz, color=10, linestyle = 2

plots, [2e4,1e5],[190,190]
plots, [2e4,1e5],[180,180], linestyle = 2
xyouts, 1e5, 190, " PCE"
xyouts, 1e5, 180, " +Diff"

xyouts, 2e4, 160,"NO", color = 2
xyouts, 2e4, 155, "N(4S)", color = 10

