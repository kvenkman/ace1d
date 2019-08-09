n_window = 0
scale = 1

    no_loss_t1 = zpid.i_no
    no_loss_t2 = zpid.j_no
    no_loss_t3 = reform(coeffmatrix[spindex.no,spindex.n2_p, *])*zionnow.n2_p;#
    no_loss_t4 = reform(coeffmatrix[spindex.no,spindex.n_p , *])*zionnow.n_p ; #
    no_loss_t5 = reform(coeffmatrix[spindex.no,spindex.o_p , *])*zionnow.o_p ; #
    no_loss_t6 = reform(coeffmatrix[spindex.no,spindex.o2_p, *])*zionnow.o2_p  ; # 
    no_loss_t7 = reform(coeffmatrix[spindex.no,spindex.n2p,  *])*zminornow.n2p ; # ; The rate coefficient is zero for this process
    no_loss_t8 = reform(coeffmatrix[spindex.no,spindex.n2d , *])*zminornow.n2d ; # 
    no_loss_t9 = reform(coeffmatrix[spindex.no,spindex.n4s , *])*zminornow.n4s ; #
    no_loss_t10 = reform(coeffmatrix[spindex.no,spindex.o2d_p , *])*zionnow.o2d_p
    
    no_loss_total = no_loss_t1 + no_loss_t2 + no_loss_t3 + no_loss_t4 + no_loss_t5 + $
    				no_loss_t6 + no_loss_t7 + no_loss_t8 + no_loss_t9 + no_loss_t10

    no_prod_t1 = reform(ratematrix[spindex.n2p,spindex.o2, *]) ; # 
    no_prod_t2 = reform(ratematrix[spindex.n2d,spindex.o2, *]) ; #
    no_prod_t3 = reform(ratematrix[spindex.n4s,spindex.o2, *]) ; #
    no_prod_t4 = reform(ratematrix[spindex.n2a,spindex.o2, *]) ; Will be zero, unless N2(A) densities are defined. ; #
    no_prod_t5 = reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p,spindex.o2,spindex.no] ; #

    no_prod_total = no_prod_t1 + no_prod_t2 + no_prod_t3 + no_prod_t4 + no_prod_t5

window, n_window, xsize=scale*639, ysize=scale*512
    
cgplot, no_loss_t1, model.zp, /xl, xr = [1e-12, 1]
cgoplot, no_loss_t2, model.zp, color = 'red'
cgoplot, no_loss_t3, model.zp, color = 'blue'
cgoplot, no_loss_t4, model.zp, color = 'orange'
cgoplot, no_loss_t5, model.zp, color = 'green'
cgoplot, no_loss_t6, model.zp, color = 'red', linestyle = 2
cgoplot, no_loss_t7, model.zp, color = 'blue', linestyle = 2
cgoplot, no_loss_t8, model.zp, color = 'green', linestyle = 2
cgoplot, no_loss_t9, model.zp, color = 'orange', linestyle = 2
cgoplot, no_loss_t10, model.zp, linestyle = 2

cgoplot, no_loss_total, model.zp, thick = 3

window, n_window+1, xsize=scale*639, ysize=scale*512
cgplot, no_prod_t1, model.zp, /xl, xr = [1e-10, 1e5]
cgoplot, no_prod_t2, model.zp, color = 'red'
cgoplot, no_prod_t3, model.zp, color = 'blue'
cgoplot, no_prod_t4, model.zp, color = 'orange'
cgoplot, no_prod_t5, model.zp, color = 'green'

cgoplot, no_prod_total, model.zp, thick = 3

