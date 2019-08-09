pro ACE_1D_pcechemistry,zmajnow,zminornow,zionnow,zmaj,zminor,zion,zpid,zei,model,chema,chemf,chemg,chemk,chemq,chemr, p_ave,$
ratea,rateg,ratek,rateq,rater,termsn2p,termsn2d

; solve all the PCE expressions for each species

 ;;;;;;N2P;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
P1_2P= 2. *zei.PED_N2 * zmajnow.n2den * (1.-chemf.F_1S-chemf.F_1D)
P2_2P= 2. *zion.N2_P  *zion.ELEC*chema.A_2*chemf.F_18
P3_2P= 2. *zmajnow.N2den *zpid.J_N2*(1.-chemf.F_3S-chemf.F_3D)
P4_2P= 2. *zei.AED_N2 * zmajnow.n2den *chemf.F_20
P5_2P= .25*zmajnow.N2den *zpid.DI_N2
P6_2P= zminor.N_4S*0.;pee_n2p
P7_2p = zei.di_n2*zmajnow.n2den
; old p7_20
;fix120=fltarr(zmaj.tpts);these two lines account for the relative independence of the scaling (i.e. constant factor of 6) with altitude below 120 km.  
;fix120(40:*)=1.
;P7_2P=.25*zei.PEI_N2/[6.+25.*((zmaj.alt/1.e5-120)/120.)*fix120];this is the pedin2
;P7_2P=smooth(P7_2P,2)


L1_2P= ace_ratio(rateg.g_91 , zminornow.n_2p)    ; zmaj.O2den*chemg.G_91
L2_2P= ace_ratio(rateg.g_90 , zminornow.n_2p)    ; zmaj.Oden*chemg.G_90
L3_2P= ace_ratio(rateq.q_93 , zminornow.n_2p)    ; zminor.NO*chemq.q_93
L4_2P= ace_ratio(rateq.q_92 , zminornow.n_2p)    ; zmaj.N2den*chemq.q_92
L5_2P= ace_ratio(rateq.q_96 , zminornow.n_2p)    ; zion.elec*chemq.q_96
L6_2P= rater.R_95 + rater.R_94

PTOT_N2P=P1_2P+P2_2P+P3_2P+P4_2P+P5_2P+P6_2P+P7_2P
LTOT_N2P=L1_2P+L2_2P+L3_2P+L4_2P+L5_2P+L6_2P
 
termsn2p={p1:p1_2p,p2:p2_2p,p3:p3_2p,p4:p4_2p,p5:p5_2p,p6:p6_2p,p7:p7_2p,ptot:ptot_n2p,$
l1:l1_2p,l2:l2_2p,l3:l3_2p,l4:l4_2p,l5:l5_2p,l6:l6_2p,ltot:ltot_n2p}
  
;zminor.n_2p=termsn2p.ptot/termsn2p.ltot
print,'PCE N(2P)'
zminor.n_2p = ace_1d_pcecalc(zminornow.n_2p,termsn2p.ptot,termsn2p.ltot,model.del_time)


;******************************** N_2D************************************************
P1_2D= 2.*zei.PED_N2*zmajnow.n2den * chemf.F_1D;*F_16
P2_2D= zion.N2_P*zmajnow.Oden*chemk.K_20*chemf.F_4
P3_2D= 2.*zion.N2_P*zion.ELEC*chema.A_2*chemf.F_12
P4_2D= 2.*zmajnow.N2den*zpid.J_N2*chemf.F_3d
P5_2D= zion.NO_P*zion.ELEC*chemf.F_2*chema.A_3
P6_2D= 2.*zei.AED_N2*zmajnow.n2den*chemf.F_1d*chemf.F_20;replace f_1 with f_20*f_21 when the n2p is running because the auroral pe spec is different than the dayglow

;**********NO TEMP DEPENDENCE***************  
P7_2D=zei.PEI_n2*zmajnow.n2den*.57*zmajnow.oden*[[[[1./[[p_ave(0,3)+p_ave(1,3)*exp((-1.*(zmajnow.zz-p_ave(2,3))^2)/(2.*p_ave(3,3)^2))]]]/[chemq.q_80(3)*zmajnow.oden+chemq.q_81(3)*zmajnow.o2den+chema.A_8(3)]]*chemq.q_80(3)]+$
[[[1./[[p_ave(0,4)+p_ave(1,4)*exp((-1.*(zmajnow.zz-p_ave(2,4))^2)/(2.*p_ave(3,4)^2))]]]/[chemq.q_80(4)*zmajnow.oden+chemq.q_81(4)*zmajnow.o2den+chema.A_8(4)]]*chemq.q_80(4)]+$
[[[1./[[p_ave(0,5)+p_ave(1,5)*exp((-1.*(zmajnow.zz-p_ave(2,5))^2)/(2.*p_ave(3,5)^2))]]]/[chemq.q_80(5)*zmajnow.oden+chemq.q_81(5)*zmajnow.o2den+chema.A_8(5)]]*chemq.q_80(5)]+$
[[[1./[[p_ave(0,6)+p_ave(1,6)*exp((-1.*(zmajnow.zz-p_ave(2,6))^2)/(2.*p_ave(3,6)^2))]]]/[chemq.q_80(6)*zmajnow.oden+chemq.q_81(6)*zmajnow.o2den+chema.A_8(6)]]*chemq.q_80(6)]];+$

P7_2D=p7_2d+zei.aEI_n2*zmajnow.n2den*.57*zmajnow.oden*[[[[1./[[p_ave(0,3)+p_ave(1,3)*exp((-1.*(zmajnow.zz-p_ave(2,3))^2)/(2.*p_ave(3,3)^2))]]]/[chemq.q_80(3)*zmajnow.o+chemq.q_81(3)*zmajnow.o2+chema.A_8(3)]]*chemq.q_80(3)]+$
[[[1./[[p_ave(0,4)+p_ave(1,4)*exp((-1.*(zmajnow.zz-p_ave(2,4))^2)/(2.*p_ave(3,4)^2))]]]/[chemq.q_80(4)*zmajnow.oden+chemq.q_81(4)*zmajnow.o2den+chema.A_8(4)]]*chemq.q_80(4)]+$
[[[1./[[p_ave(0,5)+p_ave(1,5)*exp((-1.*(zmajnow.zz-p_ave(2,5))^2)/(2.*p_ave(3,5)^2))]]]/[chemq.q_80(5)*zmajnow.oden+chemq.q_81(5)*zmajnow.o2den+chema.A_8(5)]]*chemq.q_80(5)]+$
[[[1./[[p_ave(0,6)+p_ave(1,6)*exp((-1.*(zmajnow.zz-p_ave(2,6))^2)/(2.*p_ave(3,6)^2))]]]/[chemq.q_80(6)*zmajnow.oden+chemq.q_81(6)*zmajnow.o2den+chema.A_8(6)]]*chemq.q_80(6)]];+$

;P7_2D=P6_2D*0.
P8_2D= zion.N_P*zmajnow.O2den*chemk.K_41*chemf.F_5*chemf.F_13
P9_2D= rateg.g_90 * chemf.f_17    ; zminor.N_2P*zmajnow.Oden*chemg.G_90*chemf.F_17*.47;with.03 to n4s=.5-the other half to NO+
P10_2D=.25*zmajnow.N2den*zpid.DI_N2;siskind 1995 assume .3 doublet yield.
P11_2D=rater.r_95   ;             zminor.N_2P*chemr.R_95
P12_2D=zminornow.N_4S*0.;pee_n2d

;P13_2D=.25*zei.PEI_N2/[6.+25.*((zmajnow.zz-120.)/120.)*fix120] ; need to check this!!!!!!!!!!!!!!
;P13_2D=smooth(P13_2D,2)
p13_2D=zmajnow.n2den*zei.di_n2 ; need to make sure this is right, there should be some branching ratio.....!!!!!

L1_2D= ace_ratio(rateg.g_15,zminornow.n_2d)   ; zmajnow.O2den*chemg.G_15
L2_2D= ace_ratio(rateq.q_56,zminornow.n_2d)   ; zion.ELEC*chemq.Q_56
L3_2D= ace_ratio(rateq.q_50,zminornow.n_2d)   ; zmajnow.Oden*chemq.Q_50
L4_2D= ace_ratio(rateg.g_35,zminornow.n_2d)   ; zminor.NO*chemg.G_35
L5_2D= rater.r_5                           ; fltarr(model.nlev) + chemr.R_5
L6_2D= ace_ratio(rateq.q_52,zminornow.n_2d)   ; zmajnow.N2den*chemq.Q_52
       
PTOT_N2D=P1_2D+P2_2D+P3_2D+P4_2D+P5_2D+P6_2D+P7_2D+P8_2D+P9_2D+P10_2D+P11_2D+P12_2D+P13_2D
LTOT_N2D=L1_2D+L2_2D+L3_2D+L4_2D+L5_2D+L6_2D

termsn2d={p1:p1_2d,p2:p2_2d,p3:p3_2d,p4:p4_2d,p5:p5_2d,p6:p6_2d,p7:p7_2d,p8:p8_2d,p9:p9_2d,p10:p10_2d,p11:p11_2d,$
p12:p12_2d,p13:p13_2d,ptot:ptot_n2d,l1:l1_2d,l2:l2_2d,l3:l3_2d,l4:l4_2d,l5:l5_2d,l6:l6_2d,ltot:ltot_n2d}
  
;zminor.N_2D=(zminor.N_2D + termsn2d.PTOT*model.del_time)/(1 + termsn2d.LTOT*model.del_time)
print,'PCE N(2D)'
zminor.n_2d = ace_1d_pcecalc(zminornow.n_2d,termsn2d.ptot,termsn2d.ltot,model.del_time)    
     
;******************** IONS****************************************

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; O+ (4S)
O_P_PR = $
zmajnow.oden*(zpid.i_o_op4s + zpid.i_o_op2p)                     + $ ; zmajnow.Oden*zpid.I_O*(1.-chemf.F_O)
zmajnow.O2den*zpid.DI_O2                                         + $ ; dissociative ionization gives 100% yield to O+(4S)
zmajnow.oden*(zei.pei_o_op4s + zei.pei_o_op2p)                   + $ ; zei.PEI_O*(1.-.3)
zmajnow.oden*zei.AEI_O                                           + $ ; it's probably not 100% yield, but this is the best we can do
ratek.K_40                                                    + $ ; zion.N_P*zmajnow.Oden*chemk.K_40
ratek.K_41*chemf.F_7                                          + $ ; zion.N_P*zmajnow.O2den*chemk.K_41*chemf.F_7
ratek.K_43*chemf.F_10                                         + $ ; zion.N_P*zminor.NO*chemk.K_43*chemf.F_10
ratek.K_20*(1.-chemf.F_4)*(1.-chemf.F_14)                     + $ ; zion.N2_P*zmajnow.Oden*chemk.K_20*(1.-chemf.F_4)*(1.-chemf.F_14) - the f values are correct, but why not have one for this channel only?
ratek.k_00
 
O_P_LR = $
ace_ratio(ratek.K_02,zion.o_p)   + $  ; zmajnow.N2den*chemk.K_02
ace_ratio(ratek.K_01,zion.o_p)   + $  ; zmajnow.O2den*chemk.K_01
ace_ratio(ratek.K_03,zion.o_p)        ; zminor.NO*chemk.K_03 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; O+2D
O2D_P_PR = $
zmajnow.oden*zpid.i_o_op2d    + $ ; zmajnow.Oden*zpid.I_O*chemf.F_O
zmajnow.oden*zei.pei_o_op2d       ; zei.pei_O*.3

O2D_P_LR = $
ace_ratio(ratek.K_72,zion.o2d_p)          + $ ; zmajnow.N2den*chemk.K_72
ace_ratio(ratek.K_71,zion.o2d_p)          + $ ; zmajnow.O2den*chemk.K_71
ace_ratio(ratek.k_00,zion.o2d_p)              ; zmajnow.Oden*5.e-10 (Is this off by a factor of 10?)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; N+ 

N_P_PR = $
zminor.N_4S*zpid.I_N4s           + $
zmajnow.N2den*zpid.DI_N2            + $
ratek.K_42                       + $ ; zion.N2_P*zminor.N_4S*chemk.K_42
zmajnow.n2den*zei.di_n2

N_P_LR = $
ace_ratio(ratek.K_41,zion.n_p)   + $ ;zmajnow.O2den*chemk.K_41
ace_ratio(ratek.K_40,zion.n_p)   + $ ; zmajnow.Oden*chemk.K_40
ace_ratio(ratek.K_43,zion.n_p)       ; zminor.NO*chemk.K_43

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; N2+
N2_P_PR = $
zmajnow.N2den*zpid.I_N2                      + $
ratek.K_72                                + $ ; zmajnow.N2den*zion.O2D_P*chemk.K_72
zmajnow.n2den*zei.PEI_N2                     + $
zmajnow.n2den*zei.AEI_N2                     + $
ratek.K_43*chemf.F_9                          ; zion.N_P*zminor.NO*chemk.K_43*chemf.F_9 (f_9=0.07)

N2_P_LR = $
ace_ratio(ratek.K_20,zion.n2_p)           + $ ; zmajnow.Oden*chemk.K_20
ace_ratio(ratek.K_21,zion.n2_p)           + $ ; zmajnow.O2den*chemk.K_21 
ace_ratio(ratea.A_2 ,zion.n2_p)           + $ ; zion.ELEC*chema.A_2
ace_ratio(ratek.K_23,zion.n2_p)           + $ ; zminor.NO*chemk.K_23 
ace_ratio(ratek.K_42,zion.n2_p)               ; zminor.N_4S*chemk.K_42

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; O2+
O2_P_PR = $
zmajnow.O2den*zpid.I_O2                      + $
ratek.K_21                                + $ ; zion.N2_P*zmajnow.O2den*chemk.K_21
zmajnow.o2den*zei.PEI_O2                     + $
zmajnow.o2den*zei.AEI_O2                     + $
ratek.K_01                                + $ ; zion.O_P*zmajnow.O2den*chemk.K_01
ratek.K_71                                + $ ; zion.O2D_P*zmajnow.O2den*chemk.K_71
ratek.K_41*chemf.F_5                          ; zion.N_P*zmajnow.O2den*chemk.K_41*chemf.F_5

O2_P_LR = $
ace_ratio(ratek.K_13,zion.o2_p)     + $ ; zminor.NO*chemk.K_13
ace_ratio(ratek.K_14,zion.o2_p)     + $ ; zminor.N_4S*chemk.K_14
ace_ratio(ratea.A_1 ,zion.o2_p)         ; zion.ELEC*chema.A_1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NO+
NO_P_PR = $
zminor.NO*zpid.I_NO                                                 + $
ratek.K_13                                                           + $ ; ; O2+ + NO -> NO+ + O
ratek.K_43*chemf.F_8                                                 + $ ; zminor.NO*zion.N_P*chemk.K_43*chemf.F_8 
ratek.K_23                                                          + $ ; zminor.NO*zion.N2_P*chemk.K_23
ratek.K_03                                                           + $ ; zminor.NO*zion.O_P*chemk.K_03
ratek.K_20*(chemf.F_4+(1.-chemf.F_4)*chemf.F_14)                     + $ ; zion.N2_P*zmajnow.Oden*chemk.K_20*(chemf.F_4+(1.-chemf.F_4)*chemf.F_14)
ratek.K_02                                                           + $
ratek.K_14                                                          + $
ratek.K_41*chemf.F_6                             + $ ; N+ + O2 -> NO+ + O   
rateg.G_90*chemf.F_17*chemf.f_170                                        ; zminor.N_2P*zmajnow.Oden*chemg.G_90*chemf.F_17*.5

NO_P_LR = $
ace_ratio(ratea.A_3,zionnow.no_p) ; zionnow.ELEC*chema.A_3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATE IONS
  ;N_2D=(N_2D + PTOT_N2D*DEL_TIME)/(1 + LTOT_N2D*DEL_TIME)
        ;N_2P=(N_2P + PTOT_N2P*DEL_TIME)/(1 + LTOT_N2P*DEL_TIME)
        
;zion.O_P     =(zionnow.O_P     + O_P_PR    *model.del_time)/(1.0 + O_P_LR     *model.del_time)
;zion.O2D_P   =(zionnow.O2D_P   + O2D_P_PR  *model.del_time)/(1.0 + O2D_P_LR   *model.del_time)
;zion.N_P     =(zionnow.N_P     + N_P_PR    *model.del_time)/(1.0 + N_P_LR     *model.del_time)
;zion.N2_P    =(zionnow.N2_P    + N2_P_PR   *model.del_time)/(1.0 + N2_P_LR    *model.del_time)
;zion.O2_P    =(zionnow.O2_P    + O2_P_PR   *model.del_time)/(1.0 + O2_P_LR    *model.del_time)
;zion.NO_P   =(zionnow.NO_P     + NO_P_PR   *model.del_time)/(1.0 + NO_P_LR    *model.del_time)

print,'PCE O+'
zion.o_p   = ace_1d_pcecalc(zionnow.o_p  , o_p_pr  , o_p_lr  , model.del_time)
print,'PCE O(2D)+'
zion.o2d_p = ace_1d_pcecalc(zionnow.o2d_p, o2d_p_pr, o2d_p_lr, model.del_time)
print,'PCE N+'
zion.n_p   = ace_1d_pcecalc(zionnow.n_p  , n_p_pr  , n_p_lr  , model.del_time)
print,'PCE N2+'
zion.n2_p  = ace_1d_pcecalc(zionnow.n2_p , n2_p_pr , n2_p_lr , model.del_time)
print,'PCE O2+'
zion.o2_p  = ace_1d_pcecalc(zionnow.o2_p , o2_p_pr , o2_p_lr , model.del_time)
print,'PCE NO+'
zion.no_p  = ace_1d_pcecalc(zionnow.no_p , no_p_pr , no_p_lr , model.del_time)

zion.ELEC = zion.O_P + zion.O2D_P + zion.N_P + zion.N2_P + zion.O2_P + zion.NO_P 

return
end