PRO ace_1d_major_chem, zmajnow, zminornow, zionnow, zei, zpid, coeffmatrix,  $
                       yieldmatrix, spindex, o_o2_terms, i, mass, model, pconst
                       
; The major solver distinguishes between processes that converts O2 to O (and vice versa)
; and other processes that may produce either species.

; Signs for the loss terms are taken care at the end
    pkt = model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn) ; Number density

; O2 loss terms that produce O
    s11 = reform($
                 coeffmatrix[spindex.n2d, spindex.o2, *]*zminornow.n2d + $
                 coeffmatrix[spindex.n2p, spindex.o2, *]*zminornow.n2p + $
                 coeffmatrix[spindex.n4s, spindex.o2, *]*zminornow.n4s + $
                 coeffmatrix[spindex.o_p, spindex.o2, *]*zionnow.o_p + $
                 coeffmatrix[spindex.o2p_p, spindex.o2, *]*zionnow.o2p_p + $
                 coeffmatrix[spindex.o2d_p, spindex.o2, *]*zionnow.o2d_p + $
                 coeffmatrix[spindex.n_p, spindex.o2, *]*zionnow.n_p*yieldmatrix[spindex.n_p,spindex.o2,spindex.o] + $
                 zpid.j_o2 + $ ; O2 + hv -> O(3P) + O(1D)
                 zpid.di_o2 + $
                 zei.di_o2 + $
                 zei.ped_o2 $
                )

; Production of O2 from O
    s12 = reform($
                 coeffmatrix[spindex.o, spindex.o, *]*zmajnow.oden*(zmajnow.o2den+zmajnow.n2den)*mass.o2/mass.o $
                )

; Production of O from O2
    s21 = (reform($
                 coeffmatrix[spindex.n2d, spindex.o2, *]*zminornow.n2d + $
                 coeffmatrix[spindex.n2p, spindex.o2, *]*zminornow.n2p + $
                 coeffmatrix[spindex.n4s, spindex.o2, *]*zminornow.n4s + $
                 coeffmatrix[spindex.o_p, spindex.o2, *]*zionnow.o_p + $
                 coeffmatrix[spindex.o2p_p, spindex.o2, *]*zionnow.o2p_p + $
                 coeffmatrix[spindex.o2d_p, spindex.o2, *]*zionnow.o2d_p + $
                 coeffmatrix[spindex.n_p, spindex.o2, *]*zionnow.n_p*yieldmatrix[spindex.n_p,spindex.o2,spindex.o]) + $
                 zpid.j_o2 + $   ; O2 + hv -> O(3P) + O(1D)
                 zpid.di_o2 + $
                 zei.di_o2 + $
                 2.*zei.ped_o2) $
                 *mass.o/mass.o2

; Loss of O producing O2
    s22 =  reform($
                 2.*coeffmatrix[spindex.o, spindex.o, *]*((zmajnow.o2den+zmajnow.n2den)*zmajnow.oden) $
                 )

; Explicit production/loss terms for O and O2:
; For O2
    s10 = reform($
                ; Production Terms
                 coeffmatrix[spindex.o2_p, spindex.no, *]*zionnow.o2_p*zminornow.no + $
                 coeffmatrix[spindex.o2_p, spindex.n2d, *]*zionnow.o2_p*zminornow.n2d*yieldmatrix[spindex.o2_p,spindex.n2d,spindex.o2] - $
                ; Loss Terms
                 coeffmatrix[spindex.n2_p, spindex.o2, *]*zionnow.n2_p*zmajnow.o2den - $
                 coeffmatrix[spindex.n_p, spindex.o2, *]*zionnow.n_p*zmajnow.o2den*(1. - yieldmatrix[spindex.n_p,spindex.o2,spindex.o]) - $
                 (zpid.i_o2 + zei.pei_o2 + zei.aei_o2)*zmajnow.o2den $
                 )*mass.o2/(zmajnow.barm*pkt)

; For O
    s20 = reform($
                 ; Production Terms
                 coeffmatrix[spindex.n4s, spindex.no, *]*zminornow.n4s*zminornow.no + $
                 coeffmatrix[spindex.o2_p, spindex.n2d, *]*zionnow.o2_p*zminornow.n2d*yieldmatrix[spindex.o2_p,spindex.n2d,spindex.o] + $
                 coeffmatrix[spindex.n2d, spindex.no, *]*zminornow.n2d*zminornow.no + $
                 coeffmatrix[spindex.n2d, spindex.o_p, *]*zminornow.n2d*zionnow.o_p + $
                 coeffmatrix[spindex.n_p, spindex.no, *]*zionnow.n_p*zminornow.no*yieldmatrix[spindex.n_p, spindex.no, spindex.o] + $
                 zpid.j_no*zminornow.no + $
                 coeffmatrix[spindex.o2_p, spindex.n4s, *]*zionnow.o2_p*zminornow.n4s + $
                 coeffmatrix[spindex.o_p, spindex.n4s, *]*zionnow.o_p*zminornow.n4s   + $
                 coeffmatrix[spindex.o_p, spindex.no, *]*zionnow.o_p*zminornow.no     + $
                 coeffmatrix[spindex.o2p_p, spindex.n2, *]*zionnow.o2p_p*zmajnow.n2den   + $
                 coeffmatrix[spindex.o2d_p, spindex.n4s, *]*zionnow.o2d_p*zminornow.n4s  + $
                 coeffmatrix[spindex.o2d_p, spindex.no, *]*zionnow.o2d_p*zminornow.no  + $
                 coeffmatrix[spindex.o2d_p, spindex.n2, *]*zionnow.o2d_p*zmajnow.n2den   + $
                 coeffmatrix[spindex.no_p, spindex.e, *]*zionnow.no_p*zionnow.e + $
                 coeffmatrix[spindex.o1d, spindex.o2, *]*zminornow.o1d*zmajnow.o2den + $
                 coeffmatrix[spindex.o1d, spindex.n2, *]*zminornow.o1d*zmajnow.n2den + $
                 coeffmatrix[spindex.o1d, spindex.o, *]*zminornow.o1d*zmajnow.oden + $
                 coeffmatrix[spindex.o1d, spindex.onebody, *]*zminornow.o1d + $
                 coeffmatrix[spindex.o2_p, spindex.e, *]*zionnow.o2_p*zionnow.e*yieldmatrix[spindex.o2_p,spindex.e,spindex.o] - $
                 ; Loss terms
                 coeffmatrix[spindex.n2p, spindex.o, *]*zminornow.n2p*zmajnow.oden*$
                 yieldmatrix[spindex.n2p, spindex.o, spindex.no_p] - $ ; This loss has a yield attached to it on purpose
                 
                 coeffmatrix[spindex.n_p, spindex.o, *]*zionnow.n_p*zmajnow.oden - $
                 coeffmatrix[spindex.n2_p, spindex.o, *]*zionnow.n2_p*zmajnow.oden - $
                 (zpid.i_o_op4s + zpid.i_o_op2p + zpid.i_o_op2d)*zmajnow.oden - $
                 (zei.pei_o_op4s + zei.pei_o_op2p + zei.pei_o_op2d)*zmajnow.oden - $
                 zei.aei_o*zmajnow.oden $
                 )*mass.o/(zmajnow.barm*pkt)
                 

    o_o2_terms.s11 = (-1.)*s11
    o_o2_terms.s12 = s12
    o_o2_terms.s21 = s21
    o_o2_terms.s22 = (-1.)*s22
    o_o2_terms.s10 = s10
    o_o2_terms.s20 = s20

END
