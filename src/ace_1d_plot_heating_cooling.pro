pro ace_1d_plot_heating_cooling,zmaj,pconst,heatterms,coolterms

window,10,xsize=1500,ysize=750
!p.multi=[0,2,1]

plot,alog10(heatterms.q_total),zmaj.zz,/nodata,xr=[-4,4],xtitle='Log!d10!n Cooling (K per Day)',ytitle='Altitude (Km)',yr=[100,600]

oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,coolterms.noc      )),zmaj.zz,color=!ct.BLUE
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,coolterms.nov      )),zmaj.zz,color=!ct.Navy
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,coolterms.co2_cool )),zmaj.zz,color=!ct.red
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,coolterms.o3p_cool )),zmaj.zz,color=!ct.green
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,coolterms.net_cool )),zmaj.zz,color=!ct.black

n=5
items=['NOC','NOV','CO2','O(3P)','Total']
clrs =[!ct.blue,!ct.navy,!ct.red,!ct.green,!ct.black]
syms=intarr(n)
ls=intarr(n)
legend,items,colors=(clrs),textcolors=(clrs),psym=syms,linestyle=ls,/left,/bottom,linsize=1.0,pspacing=0.5,box=0

plot,alog10(heatterms.q_total),zmaj.zz,/nodata,xr=[-4,4],xtitle='Log!d10!n Heating (K per Day)',ytitle='Altitude (Km)',yr=[100,600]

oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_chem))     ,zmaj.zz,color=!ct.navy
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_neutneut)) ,zmaj.zz,color=!ct.navy,linestyle=1
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_ionrec))   ,zmaj.zz,color=!ct.navy,linestyle=2
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_ionneut))  ,zmaj.zz,color=!ct.navy,linestyle=3
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_euv))      ,zmaj.zz,color=!ct.violet
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_thermale)) ,zmaj.zz,color=!ct.red
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_srb))      ,zmaj.zz,color=!ct.green
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_src))      ,zmaj.zz,color=!ct.amber
oplot,alog10(ace_1d_ergs2kperday(pconst,zmaj,heatterms.q_total))    ,zmaj.zz,color=!ct.black

n=9
items=['Net Chem','Chem: nn','Chem: ir','Chem: in','EUV','e','SRB','SRC','Total']
clrs =[!ct.navy,!ct.navy,!ct.navy,!ct.navy,!ct.violet,!ct.red,!ct.green,!ct.amber,!ct.black]
syms=intarr(n)
ls=[0,1,2,3,0,0,0,0,0]
legend,items,colors=(clrs),textcolors=(clrs),psym=syms,linestyle=ls,/left,/bottom,linsize=1.0,pspacing=0.5,box=0



!p.multi=0
return
end

