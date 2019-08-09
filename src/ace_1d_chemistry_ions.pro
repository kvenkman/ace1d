pro ace_1d_chemistry_ions,zmajnow,zminornow,zionnow,zion,$
                           ratematrix, coeffmatrix, yieldmatrix, $
                           spindex, zpid,zei,model,pconst,mass,count

; Solving for ion species

; O+ (4S)

o_p_prod = $
            zmajnow.oden*zpid.i_o_op4s + $ ;  + zpid.i_o_op2p
            zmajnow.o2den*zpid.di_o2*(0.54+0.46) + $ ; Need to check yield for these processes (D.I of O2 yielding O+(4S, 2D, 2P))
            zmajnow.oden*zei.pei_o_op4s + $ ; #
            zmajnow.o2den*zei.di_o2*(0.54+0.46) + $ ;
            zmajnow.oden*zei.aei_o + $
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
            reform(coeffmatrix[spindex.o_p,spindex.n2, *])*zmajnow.n2den + $ ;#
            reform(coeffmatrix[spindex.o_p,spindex.o2, *])*zmajnow.o2den + $ ;#
            reform(coeffmatrix[spindex.o_p,spindex.n2d, *])*zminornow.n2d + $ ;#
            reform(coeffmatrix[spindex.o_p,spindex.no, *])*zminornow.no; #

; O+(2D)

o2d_p_prod = $
            zmajnow.oden*zpid.i_o_op2d    + $ ; #
            zmajnow.oden*zei.pei_o_op2d   + $ ; #
            zmajnow.o2den*zpid.di_o2*0.24*0. + $
            zmajnow.o2den*zei.di_o2*0.24*0.  + $
            reform(ratematrix[spindex.o2p_p,spindex.e , *])*yieldmatrix[spindex.o2p_p,spindex.e , spindex.o2d_p] + $
            reform(ratematrix[spindex.o2p_p,spindex.onebody,*])*yieldmatrix[spindex.o2p_p,spindex.onebody,spindex.o2d_p] + $
            reform(ratematrix[spindex.o2p_p, spindex.n2_p,*])
                      

o2d_p_loss = $
            reform(coeffmatrix[spindex.o2d_p,spindex.n2, *])*zmajnow.n2den + $ ; #
            reform(coeffmatrix[spindex.o2d_p,spindex.o2, *])*zmajnow.o2den + $ ; #
            reform(coeffmatrix[spindex.o2d_p,spindex.o , *])*zmajnow.oden + $; #
            reform(coeffmatrix[spindex.o2d_p,spindex.e , *])*zionnow.e + $; #
            reform(coeffmatrix[spindex.o2d_p,spindex.n4s, *])*zminornow.n4s + $
            reform(coeffmatrix[spindex.o2d_p,spindex.no, *])*zminornow.no ; THIS WAS MISSING BEFORE

; O+(2P)

o2p_p_prod = $
            zmajnow.oden*zpid.i_o_op2p    + $ ; #
            zmajnow.oden*zei.pei_o_op2p   + $    ; #
            zmajnow.o2den*zpid.di_o2*0.22*0. + $
            zmajnow.o2den*zei.di_o2*0.22*0.


o2p_p_loss = $
            reform(coeffmatrix[spindex.o2p_p,spindex.n2, *])*zmajnow.n2den + $ ; #
            reform(coeffmatrix[spindex.o2p_p,spindex.o2, *])*zmajnow.o2den + $ ; #
            reform(coeffmatrix[spindex.o2p_p,spindex.o , *])*zmajnow.oden  + $; #
            reform(coeffmatrix[spindex.o2p_p,spindex.e , *])*zionnow.e     + $
            reform(coeffmatrix[spindex.o2p_p,spindex.onebody,*])


; N+ 

n_p_prod = $
            zminornow.n4s*zpid.i_n4s            + $
            zmajnow.n2den*zpid.di_n2            + $
            zmajnow.n2den*zei.di_n2             + $
            reform(ratematrix[spindex.n2_p,spindex.n4s, *])        + $ ;#
            reform(ratematrix[spindex.o2d_p,spindex.n4s , *])     + $
            reform(ratematrix[spindex.o2_p,spindex.n2d , *])*yieldmatrix[spindex.o2_p,spindex.n2d,spindex.n_p]    + $
            reform(ratematrix[spindex.o_p,spindex.n2d , *])

n_p_loss = $
            reform(coeffmatrix[spindex.n_p,spindex.o2, *])*zmajnow.o2den + $;#  
            reform(coeffmatrix[spindex.n_p,spindex.o , *])*zmajnow.oden  + $ ; #
            reform(coeffmatrix[spindex.n_p,spindex.no, *])*zminornow.no; #

; N2+

n2_p_prod = $
            zmajnow.n2den*zpid.i_n2               + $
            zmajnow.n2den*zei.pei_n2              + $
            zmajnow.n2den*zei.aei_n2              + $
            reform(ratematrix[spindex.n2,spindex.o2d_p, *]) + $; #
            reform(ratematrix[spindex.n2,spindex.o2p_p, *]) + $; #
            reform(ratematrix[spindex.n_p,spindex.no, *])*yieldmatrix[spindex.n_p,spindex.no,spindex.n2_p]   ; #

n2_p_loss = $
            reform(coeffmatrix[spindex.n2_p,spindex.o , *])*zmajnow.oden  + $ ;#
            reform(coeffmatrix[spindex.n2_p,spindex.o2, *])*zmajnow.o2den + $ ;#
            reform(coeffmatrix[spindex.n2_p,spindex.e, *])*zionnow.e      + $ ;#
            reform(coeffmatrix[spindex.n2_p,spindex.no, *])*zminornow.no  + $ ;#
            reform(coeffmatrix[spindex.n2_p,spindex.n4s, *])*zminornow.n4s ;#


; O2+

o2_p_prod = $
            zmajnow.o2den*zpid.i_o2                      + $
            zmajnow.o2den*zei.pei_o2                                + $
            zmajnow.o2den*zei.aei_o2                                + $
            reform(ratematrix[spindex.n2_p,spindex.o2, *])          + $ ;#
            reform(ratematrix[spindex.o_p,spindex.o2, *])           + $ ; #
            reform(ratematrix[spindex.o2d_p,spindex.o2, *])         + $ ; #
            reform(ratematrix[spindex.o2p_p,spindex.o2, *])         + $ ; #
            reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p, spindex.o2, spindex.o2_p] ; #



o2_p_loss = $
            reform(coeffmatrix[spindex.o2_p, spindex.no, *])*zminornow.no + $; #
            reform(coeffmatrix[spindex.o2_p, spindex.n4s, *])*zminornow.n4s + $ ; #
            reform(coeffmatrix[spindex.o2_p, spindex.n2d, *])*zminornow.n2d + $ ; #
            reform(coeffmatrix[spindex.o2_p, spindex.e, *])*zionnow.e; #


; NO+
no_p_prod = $
            zminornow.no*zpid.i_no              + $ 
            reform(ratematrix[spindex.no,spindex.o2_p,*]) + $; #
            reform(ratematrix[spindex.no,spindex.n2_p,*]) + $ ;#
            reform(ratematrix[spindex.no,spindex.n_p, *])*yieldmatrix[spindex.n_p,spindex.no,spindex.no_p] + $; #
            reform(ratematrix[spindex.no,spindex.o_p, *]) + $; #
            reform(ratematrix[spindex.no,spindex.o2d_p, *]) + $
            reform(ratematrix[spindex.n2_p,spindex.o,*])*yieldmatrix[spindex.n2_p,spindex.o,spindex.no_p] + $ ;#
            reform(ratematrix[spindex.o_p,spindex.n2, * ])  + $ ;#
            reform(ratematrix[spindex.o2_p,spindex.n4s, *]) + $ ; #
            reform(ratematrix[spindex.o2_p,spindex.n2d, *])*yieldmatrix[spindex.o2_p,spindex.n2d,spindex.no_p] + $ ; #
            reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p,spindex.o2,spindex.no_p] + $
            reform(ratematrix[spindex.n2p,spindex.o, *])*yieldmatrix[spindex.n2p,spindex.o,spindex.no_p]
            
no_p_loss = $
            reform(coeffmatrix[spindex.no_p, spindex.e, *])*zionnow.e ; #

; PCE calculations:
; 
zion.o_p_pce = ace_1d_pcecalc(zionnow.o_p  , o_p_prod  , o_p_loss  , model.del_time)
zion.o2d_p   = ace_1d_pcecalc(zionnow.o2d_p, o2d_p_prod, o2d_p_loss, model.del_time)
zion.o2p_p   = ace_1d_pcecalc(zionnow.o2p_p, o2p_p_prod, o2p_p_loss, model.del_time)
zion.n_p_pce = ace_1d_pcecalc(zionnow.n_p  , n_p_prod  , n_p_loss  , model.del_time)
zion.n2_p    = ace_1d_pcecalc(zionnow.n2_p , n2_p_prod , n2_p_loss , model.del_time);fltarr(model.nlev);
zion.o2_p    = ace_1d_pcecalc(zionnow.o2_p , o2_p_prod , o2_p_loss , model.del_time);fltarr(model.nlev);
zion.no_p    = ace_1d_pcecalc(zionnow.no_p , no_p_prod , no_p_loss , model.del_time);fltarr(model.nlev);

; Species that can undergo diffusion :
     zion.n_p = ace_1d_nplus(zionnow.n_p, n_p_prod, n_p_loss, zmajnow, zionnow, model, pconst)
      ;zion.n_p = ace_1d_pcecalc(zionnow.n_p , n_p_prod , n_p_loss , model.del_time)
      ;zion.n_p = fltarr(model.nlev)

     zion.o_p = ace_1d_oplus(zionnow.o_p, o_p_prod, o_p_loss, zmajnow, zionnow, model, pconst)
     ;zion.o_p = ace_1d_pcecalc(zionnow.o_p , o_p_prod , o_p_loss , model.del_time)
     ;zion.o_p = fltarr(model.nlev); zion.o_p_pce; 

; Assuming net charge neutrality
	zion.e = (zion.o_p + zion.o2d_p + zion.o2p_p + zion.n_p + zion.n2_p + zion.o2_p + zion.no_p)
	
END
