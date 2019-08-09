PRO ace_1d_euvac, solspec, f107d, f107a, zpid, zcol, inputs

    p = (f107d + f107a)/2.
    
    IF (P LT 80) THEN BEGIN
        solspec.current = 0.8*solspec.reference * inputs.mask
    ENDIF ELSE BEGIN
        solspec.current = solspec.reference*(1. + solspec.afac*(P-80.)) * inputs.mask
    ENDELSE
    
    solspec.current = 0.5*solspec.current ; In a global mean sense, only half of the Earth is sunlit
    
; This has been moved to solarphotonproc
;if sza lt 90.0 then begin
;  LYMAN=3.32E11*(1. + 0.18*(f107d-65.)/100.) ;SOLAR Ly-alpha
;  SIGI_NO=2.0E-18         ;LY-A X-SECTION FOR NO
;  SIG_LY_O2=1E-20         ;LY-A X-SECTION FOR O2
;  zpid.I_NO=(LYMAN*SIGI_NO*EXP(-1.*SIG_LY_O2*zcol.o2))
;endif else zpid.i_no=zpid.i_no*0.0

END
