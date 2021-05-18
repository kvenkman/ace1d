; Establish densities of ion species assuming photochemical equilibrium
; Ion species tracked in the model are
; O+/O+(4S) [Ground state]
; N+
; NO+
; O+(2P, 2D)
; O2+
; e [electron densities]

; Setting electron and ion temperatures to be the same as the neutral temperature
zion.te = zmaj.tn
zion.ti = zmaj.tn

; O+/O+(4S) densities
o_p_prod = $
            zmaj.oden*zpid.i_o_op4s + $ ;  + zpid.i_o_op2p
            zmaj.o2den*zpid.di_o2*(0.54+0.46) + $ ; Need to check yield for these processes (D.I of O2 yielding O+(4S, 2D, 2P))
            zmaj.oden*zei.pei_o_op4s + $ ; #
            zmaj.o2den*zei.di_o2*(0.54+0.46) + $ ;
            zmaj.oden*zei.aei_o + $
            reform(ratematrix[spindex.n_p,spindex.o, *]) + $ ; #
            reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p,spindex.o2,spindex.o_p] + $;#  
            reform(ratematrix[spindex.n_p,spindex.no, *])*yieldmatrix[spindex.n_p,spindex.no,spindex.o_p] + $; #
            reform(ratematrix[spindex.n2_p,spindex.o, *])*yieldmatrix[spindex.n2_p,spindex.o,spindex.o_p] + $ ;#
            reform(ratematrix[spindex.o2d_p,spindex.o, *]) + $ ;#
            reform(ratematrix[spindex.o2d_p,spindex.e, *]) + $
            reform(ratematrix[spindex.o2p_p,spindex.o, *]) + $;#
            reform(ratematrix[spindex.o2p_p,spindex.o2, *])*yieldmatrix[spindex.o2p_p,spindex.o2,spindex.o_p] + $;#
            reform(ratematrix[spindex.o2p_p,spindex.e, *])*yieldmatrix[spindex.o2p_p,spindex.e,spindex.o_p] + $
            reform(ratematrix[spindex.o2p_p,spindex.onebody,*])*yieldmatrix[spindex.o2p_p,spindex.onebody,spindex.o_p]

o_p_loss = $
            reform(coeffmatrix[spindex.o_p,spindex.n2, *])*zmaj.n2den + $ ;#
            reform(coeffmatrix[spindex.o_p,spindex.o2, *])*zmaj.o2den + $ ;#
            reform(coeffmatrix[spindex.o_p,spindex.n2d, *])*zminor.n2d + $ ;#
            reform(coeffmatrix[spindex.o_p,spindex.no, *])*zminor.no; #
STOP




; O(1D) assumed to be in PCE
p1_o1d = zpid.j_o2*zmaj.o2den; zpid.j_o2_o1d*zmaj.o2den ; The O(1D) yield for this is 1, [Fennelly and Torr, 1994]
p2_o1d = reform(ratematrix[spindex.o2_p, spindex.e, *])*yieldmatrix[spindex.o2_p,spindex.e,spindex.o1d] ;
p3_o1d = reform(ratematrix[spindex.o2, spindex.n2d, *])*yieldmatrix[spindex.n2d,spindex.o2,spindex.o1d] ; Miller and Hunter [2004] showed that this isn't a viable source of O(1D)

l1_o1d = reform(coeffmatrix[spindex.o1d, spindex.n2, *])*zmaj.n2den
l2_o1d = reform(coeffmatrix[spindex.o1d, spindex.o2, *])*zmaj.o2den
l3_o1d = reform(coeffmatrix[spindex.o1d, spindex.o, *])*zmaj.oden
l4_o1d = reform(coeffmatrix[spindex.o1d, spindex.onebody, *])

ptot_o1d = p1_o1d + p2_o1d + p3_o1d
ltot_o1d = l1_o1d + l2_o1d + l3_o1d + l4_o1d

zminor.o1d = ace_1d_pcecalc(zminor.o1d,ptot_o1d,ltot_o1d,model.del_time)
