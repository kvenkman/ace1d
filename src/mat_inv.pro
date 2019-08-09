FUNCTION mat_inv, matrix
; Inverting 2X2 matrices is easy and shouldn't be computationally expensive

    det_mat = matrix[0, 0]*matrix[1, 1] - matrix[1, 0]*matrix[0, 1]
    IF(det_mat EQ 0) THEN BEGIN
        print, 'Matrix determinant zero encountered'
        STOP
    ENDIF

    imat = fltarr(2, 2)
        imat[*, *] =  (1./det_mat)*[[matrix[1, 1],-matrix[1, 0]],[-matrix[0, 1], matrix[0, 0]]]
    return, imat
	
END
