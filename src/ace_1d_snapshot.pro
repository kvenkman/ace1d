
pro ace_1d_snapshot,model,zmaj,zminor,azminor,heatterms,zion,zpid,zei,initial,pconst

window_all=0
window_oddn=1
window_ion=2
window_tn=3
window_converge=4
window_pce_1=6
window_pce_2=7
window_heating=8
window_cooling=9

window,window_all,xs=1800,ys=768

wid=.18
sp=(1.-3.*wid)/4.
x0=sp
x1=sp+wid
x2=sp+wid+sp
x3=x2+wid
x4=x3+sp
x5=x4+wid

nytix=7
ytix=[-6,-4,-2,0,2,4,6]
y2tix=interpol(zmaj.zz,model.zp,ytix)

pos_neut  =[x0,0.15,x1,0.85]
pos_photon=[x2,0.15,x3,0.85]
pos_pe    =[x4,0.15,x5,0.85]
; neutral density and temperature

plot,alog10(zmaj.n2den),model.zp,/nodata,xtitle='Log!d10!n Density (cm!u-3!n)',ytitle='Pressure Level',yr=[model.zbot,model.ztop],xr=[0,14],xs=8,ys=8+1,position=pos_neut,$
ytickv=ytix
oplot,alog10(zmaj.n2den),model.zp,color=!ct.navy
oplot,alog10(zmaj.o2den),model.zp,color=!ct.blue
oplot,alog10(zmaj.oden) ,model.zp,color=!ct.red
oplot,alog10(zminor.n4s),model.zp,color=!ct.amber
oplot,alog10(zminor.no),model.zp,color=!ct.violet
; ions
oplot,alog10(zion.o_p),model.zp,color=!ct.red,linestyle=1
oplot,alog10(zion.o2_p),model.zp,color=!ct.blue,linestyle=1
oplot,alog10(zion.o2d_p),model.zp,color=!ct.red,linestyle=2
oplot,alog10(zion.n2_p),model.zp,color=!ct.navy,linestyle=1
oplot,alog10(zion.n_p),model.zp,color=!ct.amber,linestyle=1
oplot,alog10(zion.no_p),model.zp,color=!ct.violet,linestyle=1
oplot,alog10(zion.e),model.zp,color=!ct.grey,linestyle=1

temp0=0
temp1=2000
x0=!x.crange[0]
x1=!x.crange[1]
temp=zmaj.tn/(temp1-temp0) *(x1-x0) + x0
oplot,temp,model.zp,color=!ct.black
axis,xaxis=1,xr=[temp0,temp1],/xs,xtitle='Temperature (K)'

axis,yaxis=1,yticks=nytix-1,ytickname=string(fix(y2tix),format='$(i3)'),/ys,ytitle='Approximate Altitude (Km)'

items=['T ','O ','O!d2!n','N!d2!n','N(!u4!nS)','NO','O!u+!n','O!d2!u+!n','O!d2!u+!n D','N!d2O!u+!n','N!u+!n','NO!u+!n','e!u-!n']
clrs =[!ct.black,!ct.red,!ct.blue,!ct.navy,!ct.amber,!ct.violet,!ct.red,!ct.blue,!ct.blue,!ct.navy,!ct.amber,!ct.violet,!ct.grey]
ls=[0,0,0,0,0,0,0,0,1,0,0,0,0]+1
legend,items,colors=clrs,textcolors=clrs,/right,/top,psym=fltarr(n_elements(items)),linestyle=ls,linsize=3.0,pspacing=0.1,box=0

; Photo Ionization and Dissociation
plot,alog10(zmaj.n2den*zpid.i_n2),model.zp,/nodata,xtitle='Log!d10!n Rate (cm!u-3!n s!u-1!n)',ytitle='Pressure Level',yr=[model.zbot,model.ztop],ys=8+1,$
position=pos_photon,/noerase,xr=[0,6],title=' Photo: Ionization & Dissociation'
oplot,alog10(zmaj.oden*zpid.i_o),model.zp,color=!ct.red
oplot,alog10(zmaj.oden*zpid.i_o_op2p),model.zp,color=!ct.red,linestyle=1
oplot,alog10(zmaj.oden*zpid.i_o_op2d),model.zp,color=!ct.red,linestyle=2
oplot,alog10(zmaj.oden*zpid.i_o_op4s),model.zp,color=!ct.red,linestyle=3
oplot,alog10(zmaj.o2den*zpid.i_o2),model.zp,color=!ct.blue
oplot,alog10(zmaj.n2den*zpid.i_n2),model.zp,color=!ct.navy
oplot,alog10(zminor.no*zpid.i_no),model.zp,color=!ct.violet
oplot,alog10(zminor.n4s*zpid.i_n4s),model.zp,color=!ct.amber
oplot,alog10(zmaj.o2den*zpid.j_o2),model.zp,color=!ct.blue,psym=-1
oplot,alog10(zmaj.n2den*zpid.j_n2),model.zp,color=!ct.navy,psym=-1
oplot,alog10(zmaj.o2den*zpid.di_o2),model.zp,color=!ct.blue,linestyle=1
oplot,alog10(zmaj.n2den*zpid.di_n2),model.zp,color=!ct.navy,linestyle=1

items=['I_O','I_O O!u+!n (!u2!nP)','I_O!u+!n (!u2!nD)','I_O!u+!n (!u4!nS)','I_O!d2!n','I_N!d2!n','I_NO','I_N!u4!ns','J_O!d2!n','J_N!d2!n','DI_O!d2!n','DI_N!d2!n']
clrs =[!ct.red,!ct.red,!ct.red,!ct.red,!ct.blue,!ct.navy,!ct.violet,!ct.amber,!ct.blue,!ct.navy,!ct.blue,!ct.navy]
syms=[0,0,0,0,0,0,0,0,-1,-1,0,0]
ls=[0,1,2,3,0,0,0,0,0,0,1,1]
legend,items,colors=clrs,textcolors=clrs,psym=syms,linestyle=ls,/right,/top,linsize=3.0,pspacing=0.1,box=0

axis,yaxis=1,yticks=nytix-1,ytickname=string(fix(y2tix),format='$(i3)'),/ys,ytitle='Approximate Altitude (Km)'

; Photoelectron Ionization and Dissociation
plot,alog10(zmaj.n2den*zpid.i_n2),model.zp,/nodata,xtitle='Log!d10!n Rate (cm!u-3!n s!u-1!n)',ytitle='Pressure Level',yr=[model.zbot,model.ztop],ys=8+1,$
position=pos_pe,/noerase,xr=[0,4],title='Electron: Ionization & Dissociation'
oplot,alog10(zmaj.oden*zei.pei_o),model.zp,color=!ct.red
oplot,alog10(zmaj.oden*zei.pei_o_op2p),model.zp,color=!ct.red,linestyle=1
oplot,alog10(zmaj.oden*zei.pei_o_op2d),model.zp,color=!ct.red,linestyle=2
oplot,alog10(zmaj.oden*zei.pei_o_op4s),model.zp,color=!ct.red,linestyle=3
oplot,alog10(zmaj.o2den*zei.pei_o2),model.zp,color=!ct.blue
oplot,alog10(zmaj.n2den*zei.pei_n2),model.zp,color=!ct.navy
;oplot,alog10(zminor.no*zpei.i_no),model.zp,color=!ct.violet
;oplot,alog10(zminor.n4s*zei.i_n4s),model.zp,color=!ct.black
oplot,alog10(zmaj.o2den*zei.ped_o2),model.zp,color=!ct.blue,psym=-1
oplot,alog10(zmaj.n2den*zei.ped_n2),model.zp,color=!ct.navy,psym=-1
oplot,alog10(zmaj.o2den*zei.di_o2),model.zp,color=!ct.blue,linestyle=1
oplot,alog10(zmaj.n2den*zei.di_n2),model.zp,color=!ct.navy,linestyle=1

items=['I_O','I_O O!u+!n (!u2!nP)','I_O!u+!n (!u2!nD)','I_O!u+!n (!u4!nS)','I_O!d2!n','I_N!d2!n','J_O!d2!n','J_N!d2!n','DI_O!d2!n','DI_N!d2!n']
clrs =[!ct.red,!ct.red,!ct.red,!ct.red,!ct.blue,!ct.navy,!ct.blue,!ct.navy,!ct.blue,!ct.navy]
syms=[0,0,0,0,0,0,-1,-1,0,0]
ls=[0,1,2,3,0,0,0,0,1,1]
legend,items,colors=clrs,textcolors=clrs,psym=syms,linestyle=ls,/right,/top,linsize=3.0,pspacing=0.1,box=0

axis,yaxis=1,yticks=nytix-1,ytickname=string(fix(y2tix),format='$(i3)'),/ys,ytitle='Approximate Altitude (Km)'

window,window_oddn,xsize=768,ysize=512

plot,zminor.no,zmaj.zz,/nodata,xtitle='Density (cm!u-3!n)',ytitle='Approximate Altitude (km)',title='Odd Nitrogen',yr=[80,200],/xl,xr=[1,1e9],/ys
oplot,zminor.no,zmaj.zz,color=!ct.violet
oplot,initial.no,initial.zz,linestyle=1,color=!ct.violet
oplot,zminor.n2d,zmaj.zz,color=!ct.red
oplot,zminor.n2p,zmaj.zz,color=!ct.blue
oplot,zminor.n4s,zmaj.zz,color=!ct.navy

items=['NO','N(!u4!nS)','N(!u2!nD)','N(!u2!nP)']
clrs=[!ct.violet,!ct.navy,!ct.red,!ct.blue]
legend,items,colors=clrs,textcolors=clrs,psym=fltarr(n_elements(items)),linestyle=fltarr(n_elements(items)),/right,/top,linsize=3.0,pspacing=0.1,box=0

window,window_ion,xsize=768,ysize=512

plot,zion.e,zmaj.zz,/nodata,xtitle='Density (cm!u-3!n)',ytitle='Approximate Altitude (km)',title='Ions',yr=[80,400],/xl,xr=[1e2,1e7],/ys
oplot,(zion.o_p),zmaj.zz,color=!ct.red,linestyle=0
oplot,(zion.o2_p),zmaj.zz,color=!ct.blue,linestyle=0
oplot,(zion.o2d_p),zmaj.zz,color=!ct.red,linestyle=1
oplot,(zion.n2_p),zmaj.zz,color=!ct.navy,linestyle=0
oplot,(zion.n_p),zmaj.zz,color=!ct.amber,linestyle=0
oplot,(zion.no_p),zmaj.zz,color=!ct.violet,linestyle=0
oplot,(zion.e),zmaj.zz,color=!ct.grey,linestyle=0

items=['O!u+!n','O!d2!u+!n','O!u+!n !u2!nD','N!d2!u+!n','N!u+!n','NO!u+!n','e!u-!n']
clrs =[!ct.red,!ct.blue,!ct.red,!ct.navy,!ct.amber,!ct.violet,!ct.grey]
ls=[0,0,1,0,0,0,0]
legend,items,colors=clrs,textcolors=clrs,/right,/top,psym=fltarr(n_elements(items)),linestyle=ls,linsize=3.0,pspacing=0.1,box=0

; Another ion plot, but with IRI
ace_1d_plot_ionosphere,zmaj,zion,model

; temperature
window,window_tn,xsize=768,ysize=512

cgplot, zmaj.tn, zmaj.zz, yr = [100, 500],/nodata
cgoplot, zmaj.tn, zmaj.zz,color=!ct.red
cgoplot,initial.tn,initial.zz,color=!ct.grey

; convergence
window,window_converge,xsize=768,ysize=512
cgplot,  azminor.no[25] /max( azminor.no[25]), /xs, title = 'NO',color=!ct.violet
cgoplot, azminor.n4s[25]/max(azminor.n4s[25]), /xs, title = 'N(4S)',color=!ct.navy

; PCE assumption 1
window,window_pce_1,xsize=768,ysize=512
cgplot, zminor.n4s_pce,zmaj.zz, /xl, /ys,xr = [1e4,1e10], color = 10, yr = [80, 250];, title='Km divided by 1e3'
cgoplot, zminor.no_pce,zmaj.zz, color=2
cgoplot, initial.no, zmaj.zz,color = 2, linestyle = 1
plots, [2e4,1e5],[210,210]
plots, [2e4,1e5],[180,180], linestyle = 2
xyouts, 1e5, 210, " PCE"
xyouts, 1e5, 180, " +Diff"

xyouts, 2e4, 150,"NO", color = 2
xyouts, 2e4, 115, "N(4S)", color = 10

; PCE assumption 2
window,window_pce_2,xsize=768,ysize=512
cgplot, zminor.n4s_pce,zmaj.zz, /xl, /ys,xr = [1,1e12], color = 10, linestyle = 2, yr = [80, 500]
cgoplot, zminor.n4s,zmaj.zz, color=10
cgoplot, zminor.no_pce,zmaj.zz, color=2, linestyle = 2;, /xl, /ys,xr = [1e4,1e12],title='NO densities',ytitle='Altitude (kms)',xtitle='(cm!E-3!N)', yr = [100, 200]
cgoplot, zminor.no,zmaj.zz, color=2
xyouts, 1e10,140,'NO PCE',color = 2
xyouts, 5e3,140,'N(4S) PCE',color = 10
xyouts, 1e7,190,'NO',color = 2
xyouts, 6e7,160,'N(4S)',color = 10

; heating rates
window,window_heating,xsize=768,ysize=512
mx=max(heatterms.q_src)
plot, ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_euv),zmaj.zz,xr=[0.001,200],xtitle='Heating Rate', ytitle='Altitude (km)',yr=[80,150]
oplot,ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_euv),zmaj.zz,color=!ct.black
oplot,ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_neutneut),zmaj.zz,color=!ct.navy
oplot,ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_ionneut),zmaj.zz,color=!ct.blue
oplot,ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_ionrec),zmaj.zz,color=!ct.red
oplot,ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_srb),zmaj.zz,color=!ct.green
oplot,ace_1d_ergs2kperday(pconst,zmaj, heatterms.q_src),zmaj.zz,color=!ct.lime

return
end
