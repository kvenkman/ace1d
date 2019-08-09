PRO ace_1d_yieldmatrix, r1_index, r2_index, p_index, yield, yieldmatrix

    yieldmatrix[r1_index, r2_index, p_index] = yield
    yieldmatrix[r2_index, r1_index, p_index] = yield
    
END
