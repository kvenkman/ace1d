pro ace_1d_plot_ionosphere,zmaj,zion,model

window_ioniri=5

openr,lun,'iri_noon_equator_solstice_f10_150.txt',/get_lun
readf,lun,n_text,n_iri
s='a string'
for i=0,n_text-1 do readf,lun,s
iri={zz:fltarr(n_iri),o_p:fltarr(n_iri),o2_p:fltarr(n_iri),n_p:fltarr(n_iri),no_p:fltarr(n_iri),h_p:fltarr(n_iri),he_p:fltarr(n_iri),elec:fltarr(n_iri),n_iri:n_iri}
for i=0,n_iri-1 do begin
readf,lun,a,b,c,d,e,f,g,h
iri.zz[i]=a
iri.elec[i]=b/1e6 ; convert from m^-3 to cm^-3
iri.o_p[i]=c/100.*iri.elec[i]
iri.h_p[i]=d/100.*iri.elec[i]
iri.he_p[i]=e/100.*iri.elec[i]
iri.o2_p[i]=f/100.*iri.elec[i]
iri.no_p[i]=g/100.*iri.elec[i]
iri.n_p[i]=h/100.*iri.elec[i]
endfor
free_lun,lun

window,window_ioniri,xsize=768,ysize=512

plot,zion.e,zmaj.zz,/nodata,xtitle='Density (cm!u-3!n)',ytitle='Approximate Altitude (km)',title='Ions',yr=[80,600],/xl,xr=[1,1e7],/ys
oplot,(zion.o_p),zmaj.zz,color=!ct.red,linestyle=0
oplot,(zion.o2_p),zmaj.zz,color=!ct.blue,linestyle=0
oplot,(zion.o2d_p),zmaj.zz,color=!ct.red,linestyle=1
oplot,(zion.n2_p),zmaj.zz,color=!ct.navy,linestyle=0
oplot,(zion.n_p),zmaj.zz,color=!ct.amber,linestyle=0
oplot,(zion.no_p),zmaj.zz,color=!ct.violet,linestyle=0
oplot,(zion.e),zmaj.zz,color=!ct.grey,linestyle=0

oplot,iri.o_p ,iri.zz,psym=1,color=!ct.red
oplot,iri.o2_p,iri.zz,psym=1,color=!ct.blue
oplot,iri.n_p ,iri.zz,psym=1,color=!ct.amber
oplot,iri.no_p,iri.zz,psym=1,color=!ct.violet
oplot,iri.elec,iri.zz,psym=1,color=!ct.grey

items=['O!u+!n','O!d2!u+!n','O!u+!n !u2!nD','N!d2!u+!n','N!u+!n','NO!u+!n','e!u-!n']
clrs =[!ct.red,!ct.blue,!ct.red,!ct.navy,!ct.amber,!ct.violet,!ct.grey]
ls=[0,0,1,0,0,0,0]
legend,items,colors=clrs,textcolors=clrs,/right,/top,psym=fltarr(n_elements(items)),linestyle=ls,linsize=3.0,pspacing=0.1,box=0

return
end
