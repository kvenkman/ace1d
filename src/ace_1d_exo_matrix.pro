PRO ace_1d_exo_matrix, r1_index, r2_index, exo, exomatrix

    exomatrix[r1_index, r2_index] = exo
    exomatrix[r2_index, r1_index] = exo
    
END
