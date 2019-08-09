PRO ACE_1D_NOPHOT,model,model_sun,zmaj,zcol,zpid
 
; Nitric Oxide photolysis according to Minschwaner and Siskind param
;	In this implementation, the only opacity is from O2
;			NO and O3 opacities are not considered
;

IF model_sun.sza LT 90.0 THEN BEGIN ; SZA is not too large

	I5=3.98E11
	I9=2.21E11
	I10=2.3E11
	del5=2.3
	del9=1.5
	del10=1.5
	
	; Minschwaner and Siskind [1993], Table 1
	sigo2_5=[1.117e-23,2.447e-23,7.188e-23,3.042e-22,1.748e-21,1.112e-20]
	sigo2_9=[1.35e-22,2.991e-22,7.334e-22,3.074e-21,1.689e-20,1.658e-19]
	sigo2_10=[2.968e-22,5.831E-22,2.053e-21,8.192e-21,4.802e-20,2.655e-19]
	;
	signo_5 = fltarr(6,2)
	signo_9 = fltarr(6,2)
	signo_10 = fltarr(6,2)
	wno_5 = fltarr(6,2)
	wno_9 = fltarr(6,2)
	wno_10 = fltarr(6,2)
	;
	signo_5(0,0)=[0,1.32e-18,6.35e-19,7.09e-19,2.18e-19,4.67e-19]
	signo_5(0,1)=[0,4.41e-17,4.45e-17,4.50e-17,2.94e-17,4.35e-17]
	signo_9(0,0)=[0,0,3.05e-21,5.76e-19,2.29e-18,2.21e-18]
	signo_9(0,1)=[0,0,3.2e-21,5.71e-17,9.09e-17,6.0e-17]
	signo_10(0,0)=[1.8e-18,1.5e-18,5.01e-19,7.2e-20,6.72e-20,1.49e-21]
	signo_10(0,1)=[1.4e-16,1.52e-16,7.0e-17,2.83e-17,2.73e-17,6.57e-18]
	;
	wno_5(0,0)=[0,.0512,.136,.165,.141,.045]
	wno_5(0,1)=[0,5.68e-3,.0152,.0183,.0157,5e-3]
	wno_9(0,0)=[0,0,1.93e-3,.0973,.0975,.0348]
	wno_9(0,1)=[0,0,2.14e-4,1.08e-2,1.08e-2,3.86e-3]
	wno_10(0,0)=[.045,.18,.225,.225,.18,.045]
	wno_10(0,1)=[5.0e-3,.02,.025,.025,.02,5e-3]
	;
	ZPATH =  FLTARR(model.nlev)			;HORIZTONAL PATH- NOW CALC'ED IN GEOM
	NOSLNT = FLTARR(model.nlev)		;NITRIC OXIDE SLANT COLUMN
	J5   =NOSLNT			;JNO IN SRB (5-0) BAND
	J9   =NOSLNT                     ;JNO IN SRB (9-0) BAND
	J10  =NOSLNT                     ;JNO IN SRB (10-0) BAND
	J55  =NOSLNT		;FOR CASE W/O NO SELF ABS.
	J99  =NOSLNT
	J1010=NOSLNT
	A5 = FLTARR(6)			;6 SRB XSEC INTERVALS per each of 3 SR bands
	A9 = A5
	A10 = A9
	A55 = A5 & A99 = A5 & A1010 = A5	;FOR CASE W/O NO SELF ABS.
	;
	; calc JNO, both with and w/o NO self absorption
	;
	one = fltarr(2) + 1.0		;for summing over NO weights
	FOR I=0,model.nlev-1 DO BEGIN
		;	 NOSLNT(I) = TOTAL(NO(I:*)*ZPATH(0:TPTS-I-1))*1E5  NOT USED HERE
		NOSLNT(I) = 0.			; NITRIC OXIDE OPACITY SET TO 0.0 HERE
		A5 = (WNO_5*SIGNO_5   *EXP(-NOSLNT(I)*SIGNO_5))#ONE
		A9 = (WNO_9*SIGNO_9   *EXP(-NOSLNT(I)*SIGNO_9))#ONE
		A10= (WNO_10*SIGNO_10*EXP(-NOSLNT(I)*SIGNO_10))#ONE
		A55= (WNO_5*SIGNO_5)#ONE
		A99= (WNO_9*SIGNO_9)#ONE
		A1010=(WNO_10*SIGNO_10)#ONE
		;
		j5(I) = I5*DEL5*total(a5*exp(-sigo2_5*zcol.so2(i)))
		j5(i)  =  j5(i)/(1.0 + 1e-18*zmaj.n2(i))	;quenching of (0,0) band
		j9(I)  =  I9 *DEL9 *total(a9   *exp(-sigo2_9 *zcol.so2(i)))
		j10(I) =  I10*DEL10*total(a10  *exp(-sigo2_10*zcol.so2(i)))
		j55(I) =  I5 *DEL5 *total(a55  *exp(-sigo2_5 *zcol.so2(i)))
		j99(I) =  I9 *DEL9 *total(a99  *exp(-sigo2_9 *zcol.so2(i)))
		j1010(I)= I10*DEL10*total(a1010*exp(-sigo2_10*zcol.so2(i)))
	ENDFOR

	J_NO_old  = J5 + J9 + J10
	; zpid.J_NO = (1.5*J5 + J9 + J10)	;Murray et al.,J.Chem.Phys.,101,62,1994
	
	; Using the value calculated from Murray et al. raises the NO peak above the expected 105 km.
	; I think the above calculations are needed only for calculations below 100 km [KV]
    ; zpid.j_no = 4.5E-6*(1.+0.11*(model_sun.f107d-65.)/165.)*exp(-1.E-8*zcol.so2^0.38) ; From TIE-GCM
	zpid.j_no = 4.5E-6 ; Barth [1992]

	J_NO1_old = J55 + J99 + J1010
	J_NO1 = 1.5*J55 + J99 +J1010	;email David Siskind,March 9, 2000

ENDIF ELSE BEGIN
	zpid.j_no=fltarr(model.nlev) ; no sun to photolyze
ENDELSE
return
end
