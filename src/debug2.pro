no_p_prod1 =             zminornow.no*zpid.i_no               
no_p_prod2 =             reform(ratematrix[spindex.no,spindex.o2_p,*]) ; #
no_p_prod3 =             reform(ratematrix[spindex.no,spindex.n2_p,*])  ;#
no_p_prod4 =             reform(ratematrix[spindex.no,spindex.n_p, *])*yieldmatrix[spindex.n_p,spindex.no,spindex.no_p] ; #
no_p_prod5 =             reform(ratematrix[spindex.no,spindex.o_p, *]) ; #
no_p_prod6 =             reform(ratematrix[spindex.no,spindex.o2d_p, *]) 
no_p_prod7 =             reform(ratematrix[spindex.n2_p,spindex.o,*])*yieldmatrix[spindex.n2_p,spindex.o,spindex.no_p]  ;#
no_p_prod8 =             reform(ratematrix[spindex.o_p,spindex.n2, * ])   ;#
no_p_prod9 =             reform(ratematrix[spindex.o2_p,spindex.n4s, *])  ; #
no_p_prod10 =            reform(ratematrix[spindex.o2_p,spindex.n2d, *])  ; #
no_p_prod11 =            reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p,spindex.o2,spindex.no_p] 
no_p_prod12 =            reform(ratematrix[spindex.n2p,spindex.o, *])*yieldmatrix[spindex.n2p,spindex.o,spindex.no_p]

no_p_prod = no_p_prod1 + no_p_prod2 + no_p_prod3 + no_p_prod4 + no_p_prod5 + no_p_prod6 + no_p_prod7 + no_p_prod8 + $
			no_p_prod9 + no_p_prod10 + no_p_prod11 + no_p_prod12
			
no_p_loss = $
reform(coeffmatrix[spindex.no_p, spindex.e, *])*zionnow.e ; #

			
cgplot, no_p_prod, model.zp, /xl
cgoplot, no_p_prod1, model.zp, color = 'red'
cgoplot, no_p_prod2, model.zp, color = 'blue'
cgoplot, no_p_prod3, model.zp, color = 'green'
cgoplot, no_p_prod4, model.zp, color = 'orange'
cgoplot, no_p_prod5, model.zp, color = 'yellow'
cgoplot, no_p_prod6, model.zp, color = 'purple'

cgoplot, no_p_prod7, model.zp, color = 'red', linestyle = 2
cgoplot, no_p_prod8, model.zp, color = 'blue', linestyle = 2
cgoplot, no_p_prod9, model.zp, color = 'green', linestyle = 2
cgoplot, no_p_prod10, model.zp, color = 'orange', linestyle = 2
cgoplot, no_p_prod11, model.zp, color = 'yellow', linestyle = 2
cgoplot, no_p_prod12, model.zp, color = 'purple', linestyle = 2

