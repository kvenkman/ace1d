PRO ace_1d_chemistry_ox, zmajnow, zminornow, zionnow, zminor, zpid, coeffmatrix, $
	ratematrix, yieldmatrix, branching, spindex, model

; O(1D) assumed to be in PCE
p1_o1d = zpid.j_o2*zmajnow.o2den; zpid.j_o2_o1d*zmajnow.o2den ; The O(1D) yield for this is 1, [Fennelly and Torr, 1994]
p2_o1d = reform(ratematrix[spindex.o2_p, spindex.e, *])*yieldmatrix[spindex.o2_p,spindex.e,spindex.o1d] ;
p3_o1d = reform(ratematrix[spindex.o2, spindex.n2d, *])*yieldmatrix[spindex.n2d,spindex.o2,spindex.o1d] ; Miller and Hunter [2004] showed that this isn't a viable source of O(1D)

l1_o1d = reform(coeffmatrix[spindex.o1d, spindex.n2, *])*zmajnow.n2den
l2_o1d = reform(coeffmatrix[spindex.o1d, spindex.o2, *])*zmajnow.o2den
l3_o1d = reform(coeffmatrix[spindex.o1d, spindex.o, *])*zmajnow.oden
l4_o1d = reform(coeffmatrix[spindex.o1d, spindex.onebody, *])

ptot_o1d = p1_o1d + p2_o1d + p3_o1d
ltot_o1d = l1_o1d + l2_o1d + l3_o1d + l4_o1d

zminor.o1d = ace_1d_pcecalc(zminornow.o1d,ptot_o1d,ltot_o1d,model.del_time)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
END
