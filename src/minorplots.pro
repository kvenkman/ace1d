rho = (model.p0*exp(-model.zp))/(pconst.boltz*zmajnow.tn)*zmajnow.barm/pconst.avo
 ;  Q_NEUTNEUT      FLOAT     Array[57]
 ;  Q_IONREC        FLOAT     Array[57]
 ;  Q_IONNEUT       FLOAT     Array[57]
 ;  Q_CHEM          FLOAT     Array[57]
 ;  Q_EUV           FLOAT     Array[57]
 ;  Q_THERMALE      FLOAT     Array[57]
 ;  Q_SRB           FLOAT     Array[57]
 ;  Q_SRC           FLOAT     Array[57]
 ;  Q_TOTAL         FLOAT     Array[57]

window, 0
cgplot, alog10(heatterms.q_total/rho), zmajnow.zz, title = 'Heating rates', /ys
cgoplot, alog10(heatterms.q_neutneut/rho), zmajnow.zz,color='red'
cgoplot, alog10(heatterms.q_chem/rho), zmajnow.zz,color='green'
cgoplot, alog10(heatterms.q_src/rho), zmajnow.zz,color='purple'
cgoplot, alog10(heatterms.q_srb/rho), zmajnow.zz,color='blue'
cgoplot, alog10(heatterms.q_euv/rho), zmajnow.zz,color='yellow'
window, 1
cgplot, alog10(heatterms.q_total/rho), zmajnow.zz, title = 'Heating rates', /ys
cgoplot, alog10(heatterms.q_ionrec/rho), zmajnow.zz,color = 'red', title = 'Heating rates'
cgoplot, alog10(heatterms.q_ionneut/rho), zmajnow.zz,color = 'blue'
cgoplot, alog10(heatterms.q_thermale/rho), zmajnow.zz,color = 'purple'
window, 3, xsize = 500, ysize = 500
cgplot, alog10(zionnow.e), model.zp, /ys, xr = [2, 7], yr = [-7, 5]
cgoplot, alog10(zionnow.n_p), model.zp, linestyle = 2
END
