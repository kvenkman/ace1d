pro ace_1d_eionheating, zmajnow, zminornow, zionnow, $
    heattermse, cooltermse, heattermsi, cooltermsi, heatterms,$
    heatmatrix, spindex, pconst, model, count

; calculate electron and ion heating and cooling terms

; everything is calculated in eV/(cm^3 s) and then converted to ergs at the end
; also note that at the end the cooling of electrons by collions with neutrals is used to specify the appropriate heating term for the neutrals

; note that heating is due to photoelectrons and so the heating rate is done in ace_1d_solarphotonproc.pro
; note that as described in Smithro and Solomon 2008, the PE heating of ambient electrons is best described in two 
; terms due to solar euv above and below 55 nm.
heattermse.q_quench = (reform(heatmatrix[spindex.n2d, spindex.e, *]))*pconst.ev2erg
heattermse.q_total  = heattermse.q_pe0_55 + heattermse.q_pe55_105 + heattermse.q_quench

; cooling terms are taken from the Schunk and Nagy text (second edition page 2009); see pages 277-287
; note that we are including they key Earth terms, but more could be done with H2, CO, CH4, etc.
; see the S&N text if using this code for planetary atmospheres
; some terms as noted below are taken from Schunk and Nagy 1978
te_ti_diff = ((zionnow.te - zionnow.ti) > 0.)
te_tn_diff = ((zionnow.te - zmajnow.tn) > 0.)
ti_tn_diff = ((zionnow.ti - zmajnow.tn) > 0.)

tfactor = (1./sqrt(zionnow.te)) ;* ((zionnow.te - zmajnow.tn)>0.0)
; electrons on N2, rotational excitation
cooltermse.n2rot = 3.5e-14 * zionnow.e * zmajnow.n2den * tfactor

; electrons on O2, rotational excitation
cooltermse.o2rot = 5.2e-15 * zionnow.e * zmajnow.o2den * tfactor

; electrons on CO2, rotational excitation
cooltermse.co2rot = 5.8e-14 * zionnow.e * zminornow.co2 * tfactor

; electrons  on H2O, rotational excitation
a = 1.052e-8 + 6.043e-10 * alog(zmajnow.tn)
b = 4.180e-9 + 2.026e-10 * alog(zmajnow.tn)
cooltermse.h2orot = zionnow.e * zminornow.h2o * (a + b*alog(zmajnow.tn/zionnow.te)/(zionnow.te^(5./4.))) * te_tn_diff ;((zionnow.te - zmajnow.tn)>0.)

; electrons on N2, vibrational excitation
; tables for S&N pg 279
a0h=[2.025, -7.066, -8.211, -9.713, -10.353, -10.819, -10.183, -12.698, -14.710, -17.538]
b0h=[8.782e-4, 1.001e-2, 1.092e-2, 1.204e-2, 1.243e-2, 1.244e-2, 1.185e-2, 1.309e-2, 1.409e-2, 1.600e-2]
c0h=[2.954e-7, -3.066e-6, -3.369e-6, -3.732e-6, -3.850e-6, -3.771e-6, -3.570e-6, -3.952e-6, -4.249e-6, -4.916e-6]
d0h=[-9.562e-11, 4.436e-10, 4.891e-10, 5.431e-10, 5.600e-10, 5.385e-10, 5.086e-10, 5.636e-10, 6.058e-10, 7.128e-10]
f0h=[7.252e-15, -2.449e-14, -2.706e-14, -3.008e-14, -3.100e-14, -2.936e-14, -2.769e-14, -3.071e-14, -3.300e-14, -3.941e-14]
delta0h=[0.06, 0.08, 0.10, 0.10, 0.13, 0.15, 0.15, 0.15, 0.15, 0.15]

a0l=a0h-a0h & a0l[0]=-6.462
b0l=b0h-b0h & b0l[0]=3.151e-2
c0l=c0h-c0h & c0l[0]=-4.075e-5
d0l=d0h-d0h & d0l[0]=2.439e-8
f0l=f0h-f0h & f0l[0]=-5.479e-12
delta0l=delta0h-delta0h & delta0l[0]=0.14

a1h=[-3.413,-4.160, -5.193, -5.939, -8.261, -8.185, -10.823, -11.273]
b1h=[7.326e-3, 7.803e-3, 8.360e-3, 8.807e-3, 1.010e-2, 1.010e-2, 1.199e-2, 1.283e-2]
c1h=[-2.200e-6, -2.352e-6, -2.526e-6, -2.669e-6, -3.039e-6, -3.039e-6, -3.620e-6, -3.879e-6]
d1h=[3.128e-10, 3.352e-10, 3.606e-10, 3.806e-10, 4.318e-10, 4.318e-10, 5.159e-10, 5.534e-10]
f1h=[-1.702e-14, -1.828e-14, -1.968e-14, -2.073e-14, -2.347e-14, -2.347e-14, -2.810e-14, -3.016e-14]
delta1h=[0.11, 0.11, 0.12, 0.08, 0.10, 0.12, 0.09, 0.09]

a1l=a1h-a1h 
b1l=b1h-b1h 
c1l=c1h-c1h 
d1l=d1h-d1h 
f1l=f1h-f1h 
delta1l=delta1h-delta1h 

e1 = 3353.0 ; K
tvib=zmajnow.tn
maskh = (zionnow.te ge 1500.)
maskl = ((zionnow.te ge 300.) and (zionnow.te lt 1500.))

sum0=cooltermse.n2vib*0.0
sum1=sum0

factor0 = zionnow.e * zmajnow.n2den * (1.0 - exp(-e1/tvib))
factor1 = zionnow.e * zmajnow.n2den * (1.0 - exp(-e1/tvib))*exp(-e1/tvib)

logq0=fltarr(10,n_elements(zionnow.te))
logq1=logq0
q0=logq0
q1=logq1

for i=1,10 do begin

i1=i-1

logq0[i1,*] = (a0l[i1]*maskl + a0h[i1]*maskh) + (b0l[i1]*maskl + b0h[i1]*maskh)*zionnow.te + (c0l[i1]*maskl + c0h[i1]*maskh)*zionnow.te^2 + $
(d0l[i1]*maskl + d0h[i1]*maskh)*zionnow.te^3 + (f0l[i1]*maskl + f0h[i1]*maskh)*zionnow.te^4 - 16.0
q0[i1,*]=10.0^logq0[i1,*]

sum0 = sum0 + q0[i1,*]*(1 - exp(float(i) * e1 * (1.0/zionnow.te - 1.0/tvib))) 

if ((i ge 2) and (i le 9)) then begin 

i2=i-2

logq1[i2,*] = (a1l[i2]*maskl + a1h[i2]*maskh) + (b1l[i2]*maskl + b1h[i2]*maskh)*zionnow.te + (c1l[i2]*maskl + c1h[i2]*maskh)*zionnow.te^2 + $
(d1l[i2]*maskl + d1h[i2]*maskh)*zionnow.te^3 + (f1l[i2]*maskl + f1h[i2]*maskh)*zionnow.te^4 - 16.0
q1[i2,*]=10.0^logq1[i2,*]

sum1 = sum1 + q1[i2,*]*(1 - exp(float(i-1) * e1 * (1.0/zionnow.te - 1.0/tvib)))

endif
endfor

cooln2vib = (factor0 * sum0 + factor1 * sum1) * (zionnow.te gt zmajnow.tn)
cooltermse.n2vib = cooln2vib>0.0

; ; electrons on O2, vibrational excitation
logq=-19.9171 + 0.0267*zionnow.te - 3.9960e-5*zionnow.te^2 + 3.5187e-8*zionnow.te^3 - 1.9228e-11*zionnow.te^4 + 6.6865e-15*zionnow.te^5 -$
   1.4791e-18*zionnow.te^6 + 2.0127e-22*zionnow.te^7 - 1.5346e-26*zionnow.te^8 + 5.0148e-31*zionnow.te^9
q=10.0^logq

cooltermse.o2vib = zionnow.e * zmajnow.o2den * q * (1.0 - exp(2239. * (1.0/zionnow.te - 1.0/zmajnow.tn))) > 0.

; electrons on O, O fine structure
d = 5.0 + exp(-326.6/zmajnow.tn) + 3.0*exp(-227.7/zmajnow.tn)
s21 = 1.863e-11
s20 = 1.191e-11
s10 = 8.249e-16 * (zionnow.te^0.6) * exp(-227.7/zmajnow.tn)
tt=(1.0/zionnow.te - 1.0/zmajnow.tn)
cooltermse.ofine = zionnow.e * zmajnow.oden * (1.0/d) * ( s10*(1.0 - exp(98.9*tt )) + $
                                                          s20*(1.0 - exp(326.6*tt)) + $
                                                          s21*(1.0 - exp(227.7*tt)) )>0.0

; electrons on O, O(1D) excitation
d = 2.4e4 + 0.3*(zionnow.te - 1500.) - 1.947e-5*(zionnow.te - 1500.)*(zionnow.te - 4000.)
cooltermse.o1d = (1.57e-12 * zionnow.e * zmajnow.oden * exp(d*(zionnow.te - 3000.)/(3000.*zionnow.te)) * (exp(-22713.*te_tn_diff/(zionnow.te*zmajnow.tn)) - 1.0))>0.0

;test_o1d_cool = 1.07e-10*zionnow.e*zmajnow.oden*(zionnow.te^0.5)*exp(-2.27e4/zionnow.te)* $
;				(0.406+0.357e-4*zionnow.te-(0.333 + 0.183e-4*zionnow.te)*exp(-1.37e4/zionnow.te) - $
;				(0.456+0.174e-4*zionnow.te)*exp(-2.97e4/zionnow.te))

; These expressions are from Schunk and Nagy [1978], 43a-c
; electrons on N2, elastic
cooltermse.elasticn2 = (1.77e-19 * zionnow.e * zmajnow.n2den * (1.0 - 1.21e-4*zionnow.te)*zionnow.te) ; >0.0

; electrons on O2, elastic
cooltermse.elastico2 = (1.21e-18 * zionnow.e * zmajnow.o2den * (1.0 + 3.6e-2*sqrt(zionnow.te))*sqrt(zionnow.te)) ; >0.0

; electrons on O, elastic
cooltermse.elastico = (7.9e-19 * zionnow.e * zmajnow.oden * (1.0 + 5.7e-4*zionnow.te)*sqrt(zionnow.te)); >0.0

; electrons on ions
; This expression is from Schunk and Nagy [1978], Eq. 48
gamma=exp(0.57)
e    = 1.6d-19
zo   = 1; 8.0
zo2  = 1; 16.0
zno  = 1; 14.0
kerg = pconst.boltz; Boltzmann in ergs K^-1

n_e = zionnow.e
n_ion = zionnow.o_p + zionnow.o2_p + zionnow.no_p
ki_sq = 4*!pi*n_ion*(1.^2.)*(e^2.)/(pconst.boltz*zionnow.ti) > 1e-24
ke_sq = 4*!pi*n_e*(1.^2.)*(e^2.)/(pconst.boltz*zionnow.te) > 1e-24

; Roble & Rees [1975]; Eq. 56
 ln_gamma = alog(4*pconst.boltz*zionnow.te/(gamma^2.*(e^2.)*sqrt(ke_sq))) - ((ke_sq+ki_sq)/ki_sq)*alog(sqrt((ke_sq+ki_sq)/ke_sq))

;Pg. 104, Schunk and Nagy, Ionospheres [2009]
;ln_gamma = 15.;

cooltermse.ei = 3.2e-8*n_e*(1./zionnow.te^1.5)*ln_gamma*(zionnow.o_p + 0.5*zionnow.o2_p + 0.53*zionnow.no_p) ; Has a (te - ti) term attached

; Cooling of ions, Rees and Roble 1975, Table 3
t1 = 6.6e-14*zmajnow.n2den + 5.8e-14*zmajnow.o2den + 0.21e-14*zmajnow.oden*sqrt(zionnow.ti + zmajnow.tn)
t2 = 5.45e-14*zmajnow.o2den + 5.9e-14*zmajnow.n2den + 4.5e-14*zmajnow.oden
t3 = 5.8e-14*zmajnow.n2den + 4.4e-14*zmajnow.oden + 0.14e-14*zmajnow.o2den*sqrt(zionnow.ti + zmajnow.tn)

cool_in = (zionnow.o_p*t1 + zionnow.no_p*t2 + zionnow.o2_p*t3)*pconst.ev2erg

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; change units and sum up the cooling terms
cooltermse.n2rot  = cooltermse.n2rot  * pconst.ev2erg
cooltermse.o2rot  = cooltermse.o2rot  * pconst.ev2erg
cooltermse.co2rot = cooltermse.co2rot * pconst.ev2erg
cooltermse.h2orot = cooltermse.h2orot * pconst.ev2erg

cooltermse.n2vib  = cooltermse.n2vib  * pconst.ev2erg
cooltermse.o2vib  = cooltermse.o2vib  * pconst.ev2erg
cooltermse.ofine  = cooltermse.ofine  * pconst.ev2erg
cooltermse.o1d    = cooltermse.o1d    * pconst.ev2erg

cooltermse.elasticn2 = cooltermse.elasticn2 * pconst.ev2erg
cooltermse.elastico2 = cooltermse.elastico2 * pconst.ev2erg
cooltermse.elastico  = cooltermse.elastico  * pconst.ev2erg
cooltermse.ei        = cooltermse.ei        * pconst.ev2erg

; These don't have any temperature difference terms attached
cooltermse.explicit = cooltermse.n2vib + cooltermse.o2vib + cooltermse.ofine + cooltermse.o1d

; This has a (te - tn) term attached to it
cooltermse.neutrals_implicit = (cooltermse.n2rot +  cooltermse.o2rot +  cooltermse.co2rot +  cooltermse.h2orot + $
                      cooltermse.elasticn2 + cooltermse.elastico2 + cooltermse.elastico )
                      
; This has a (te - ti) term attached to it
cooltermse.ions_implicit = cooltermse.ei

cooltermse.total = (cooltermse.n2rot + cooltermse.o2rot + cooltermse.co2rot + cooltermse.h2orot + $
                    cooltermse.elasticn2 + cooltermse.elastico2 + cooltermse.elastico)*te_tn_diff + $
                    cooltermse.n2vib + cooltermse.o2vib + cooltermse.ofine + cooltermse.o1d + $
                    cooltermse.ei*te_ti_diff

; Cooling of electrons by collisions with ions to be heating of heattermsi
heattermsi.ei = cooltermse.ei > 1e-24

; Cooling of ions by collisions with neutrals:
cooltermsi.neut = cool_in  > 1e-24

; set cooling of electrons by neutrals to heating of neutrals by electrons
; These expressions are used only for plotting purposes. Calculations of TN treat these terms differently
heatterms.q_thermale = cooltermse.explicit + cooltermse.neutrals_implicit*te_tn_diff; (zionnow.te - zmajnow.tn)
heatterms.q_thermali = cooltermsi.neut*ti_tn_diff; (zionnow.ti - zmajnow.tn)
; heatterms.q_total = (heatterms.q_total); + heatterms.q_thermale + cooltermsi.neut)
END
