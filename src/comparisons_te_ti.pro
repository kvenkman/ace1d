restore, 'ace1d_solaractivity.sav'
restore, 'tiegcm_f107_te_ti.sav'

lev = 0.5*findgen(29) - 6.75
zp = 0.25*findgen(57) - 7


cgplot, ace1d_70.zion.te, zp
cgoplot, tiegcm_te[*, 0], lev, color = 'red'

cgoplot, ace1d_245.zion.te, zp, linestyle = 2
cgoplot, tiegcm_te[*, -1], lev, color = 'red', linestyle = 2


END
