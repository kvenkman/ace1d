PRO ace_1d_euvac, solspec, f107d, f107a, zpid, zcol, inputs

    p = (f107d + f107a)/2.
    
    IF (P LT 80) THEN BEGIN
        solspec.current = 0.8*solspec.reference * inputs.mask
    ENDIF ELSE BEGIN
        solspec.current = solspec.reference*(1. + solspec.afac*(P-80.)) * inputs.mask
    ENDELSE
    
    solspec.current = 0.5*solspec.current ; In a global mean sense, only half of the Earth is sunlit
    
END
