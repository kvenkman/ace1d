PRO ace_1d_update_crh, r1_index, r2_index, r1_den, r2_den, coeff, $
                        exo, coeffmatrix, ratematrix, heatmatrix
; This routine is used to update coefficients, rates and heating from reactions

; The value for this variable has to be updated if more species are added
nsp = 23
onebody = nsp - 4;


; 4/25/16, Note -
; When using the coefficient/rate/heat matrices in subroutines, the order of the index shouldn't matter as the matrices will be symetric
; K.V

IF ((r1_index NE onebody) AND (r2_index NE onebody)) THEN BEGIN
    ; For non one body reactions
    coeffmatrix[r1_index, r2_index, *] = coeff
    coeffmatrix[r2_index, r1_index, *] = coeffmatrix[r1_index, r2_index, *]

    ratematrix[r1_index, r2_index, *] = coeffmatrix[r1_index, r2_index, *]*r1_den*r2_den
    ratematrix[r2_index, r1_index, *] = ratematrix[r1_index, r2_index, *]

    heatmatrix[r1_index, r2_index, *] = exo*ratematrix[r1_index, r2_index, *]
    heatmatrix[r2_index, r1_index, *] = heatmatrix[r1_index, r2_index, *]
ENDIF ELSE BEGIN
    coeffmatrix[r1_index, r2_index, *] = coeff
    coeffmatrix[r2_index, r1_index, *] = coeffmatrix[r1_index, r2_index, *]

    ratematrix[r1_index, r2_index, *] = coeffmatrix[r1_index, r2_index, *]*r1_den
    ratematrix[r2_index, r1_index, *] = ratematrix[r1_index, r2_index, *]

; Exothermicities for single body reactions will be zero, though
    heatmatrix[r1_index, r2_index, *] = exo*ratematrix[r1_index, r2_index, *]
    heatmatrix[r2_index, r1_index, *] = heatmatrix[r1_index, r2_index, *]
ENDELSE
END
