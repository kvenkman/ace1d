PRO	ace_1d_chemistry_oddn, zmajnow, zminornow, zionnow, $
                           zminor, zpid, zei, $
                           ratematrix, coeffmatrix, yieldmatrix, $
                           model, spindex, termsn2p, termsn2d,$
                           termsn4s,termsno,flux,k_m,mass,pconst, model_sun, i
 
    ; The order in which species are dealt with :
    ;   NO, N(4S), N(2D), N(2P)

;;;;; Nitric Oxide

	k_n2a_o = [2.8, 3.4, 3.34, 3.48, 2.08, 4.26, 4.58, 7.4]*1e-11 ; Level specific reaction rate for v = 3-6
	n2a_no_n2d_yield = 0.57

    no_prod_t1 = reform(ratematrix[spindex.n2p,spindex.o2, *]) ; # 
    no_prod_t2 = reform(ratematrix[spindex.n2d,spindex.o2, *]) ; #
    no_prod_t3 = reform(ratematrix[spindex.n4s,spindex.o2, *]) ; #
    no_prod_t4_v3 = k_n2a_o[3]*zminornow.n2a[3,*]*zmajnow.oden*n2a_no_n2d_yield ; production of NO from N2A(v=3) by reaction with O
    no_prod_t4_v4 = k_n2a_o[4]*zminornow.n2a[4,*]*zmajnow.oden*n2a_no_n2d_yield
    no_prod_t4_v5 = k_n2a_o[5]*zminornow.n2a[5,*]*zmajnow.oden*n2a_no_n2d_yield
    no_prod_t4_v6 = k_n2a_o[6]*zminornow.n2a[6,*]*zmajnow.oden*n2a_no_n2d_yield
    no_prod_t5 = reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p,spindex.o2,spindex.no] ; #


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

    no_prod = no_prod_t1 + no_prod_t2 + no_prod_t3 + no_prod_t4_v3 $
               + no_prod_t4_v4 + no_prod_t4_v5 + no_prod_t4_v6 + no_prod_t5

    no_loss = no_loss_t1 + no_loss_t2 + no_loss_t3 + no_loss_t4 + $
              no_loss_t5 + no_loss_t6 + no_loss_t7 + no_loss_t8 + $
              no_loss_t9 + no_loss_t10
              
    ;no_loss = no_loss > 1e-7
    
    termsno = {no_loss:no_loss, no_n2d_prod:no_prod_t2, no_n4s_prod:no_prod_t3}
    
    ;zminor.no_pce = no_prod/no_loss
    zminor.no_pce = ace_1d_pcecalc(zminornow.no_pce,no_prod,no_loss,model.del_time)

    ;   no_temp = ace_1d_pcecalc(zminornow.no,no_prod,no_loss,model.del_time)
        ; no_lbc = no_temp[0] ; rhs
        ;no_lbc_flux = 0.;1e-8 ; If flux is non zero, it WILL be used in place of any other LBC
    
    ;zminor.no = ace_1d_minor_solv(zminornow.no,mass.no,zmajnow,no_prod,no_loss,k_m.no,flux.no, model, mass, pconst,i) ; zminor.no_pce ; 
    ;zminor.no_mmr = ace_1d_cm3_mmr(zminor.no, mass.no, zmajnow, model, pconst) ; 
    
    zminor.no_mmr = ace_1d_minor_no(zminornow, zmajnow, no_prod, no_loss, model, mass, pconst, model_sun)
    zminor.no = ace_1d_mmr_cm3(zminor.no_mmr, mass.no, zmajnow, model, pconst)
    
    ;cgplot, zminor.no, zmajnow.zz, yr = [100, 200], /xs;, xr = [1, 1e8];, /xl
    ;wait, 0.001
    
;;;;; N(4S)
    
    n2p_n4s_yield = 1.6e-12*(zionnow.te)^0.85/(1.6e-12*(zionnow.te)^0.85 + 9.5e-9) > 0.
    n2p_n2d_yield = 9.5e-9/(1.6e-12*(zionnow.te)^0.85 + 9.5e-9) > 0.
    
    n4s_prod_t1 = 2.*zei.ped_n2*zmajnow.n2den*yieldmatrix[spindex.n2, spindex.pe, spindex.n4s]
    n4s_prod_t2 = 2.*zei.aed_n2*zmajnow.n2den*yieldmatrix[spindex.n2, spindex.ae, spindex.n4s]
    n4s_prod_t3 = zmajnow.n2den*zpid.di_n2*0.4 ; Photoelectron dissociative ionization of N2, times N(4S) branching ratio
    n4s_prod_t4 = zmajnow.n2den*zei.di_n2*0.4 ; Dissociative ionization of N2 by EUV photons, times N(4S) branching ratio
    n4s_prod_t5 = 2.*zmajnow.n2den*zpid.j_n2*yieldmatrix[spindex.n2, spindex.ph, spindex.n4s] ; Photodissociation of N2
    
    n4s_prod_t6 = reform(ratematrix[spindex.n2_p,spindex.e, *])*yieldmatrix[spindex.n2_p, spindex.e, spindex.n4s] ; 
    n4s_prod_t7 = reform(ratematrix[spindex.n2_p,spindex.o, *])*yieldmatrix[spindex.n2_p, spindex.o, spindex.n4s] ; ;#  
    n4s_prod_t8 = reform(ratematrix[spindex.no_p,spindex.e, *])*yieldmatrix[spindex.no_p, spindex.e, spindex.n4s] ; #
    n4s_prod_t9 = reform(ratematrix[spindex.n_p,spindex.o, *])  ; #   
    n4s_prod_t10 = reform(ratematrix[spindex.n_p,spindex.o2,*])*yieldmatrix[spindex.n_p,spindex.o2,spindex.n4s]
    n4s_prod_t11 = reform(ratematrix[spindex.n_p,spindex.no , *])*yieldmatrix[spindex.n_p,spindex.no,spindex.n4s]  ; #
    n4s_prod_t12 = reform(ratematrix[spindex.o_p,spindex.n2 , *]) ;#
    n4s_prod_t13 = reform(ratematrix[spindex.n2p,spindex.e  , *])*n2p_n4s_yield ; #   
    n4s_prod_t14 = reform(ratematrix[spindex.n2p,spindex.o  , *])*yieldmatrix[spindex.n2p, spindex.o, spindex.n4s]  ; #      
    n4s_prod_t15 = reform(ratematrix[spindex.n2d,spindex.n2 , *]) ; #      
    n4s_prod_t16 = reform(ratematrix[spindex.n2d,spindex.o, *]) ; #      
    n4s_prod_t17 = reform(ratematrix[spindex.n2d,spindex.e  , *]) ; #      
    n4s_prod_t18 = reform(ratematrix[spindex.n2d,spindex.onebody, *]) ; #  
    n4s_prod_t19 = reform(ratematrix[spindex.n2p,spindex.onebody, *])*yieldmatrix[spindex.n2p, spindex.onebody, spindex.n4s] ; #  
    n4s_prod_t20 = zminornow.no*zpid.j_no ; NO photolysis    
    n4s_prod_t21 = reform(ratematrix[spindex.n2p,spindex.o2_p, *]) ; #
    
    n4s_loss_t1 = reform(coeffmatrix[spindex.n4s,spindex.n2_p, *])*zionnow.n2_p ;#
    n4s_loss_t2 = reform(coeffmatrix[spindex.n4s,spindex.o2_p, *])*zionnow.o2_p ; #    
    n4s_loss_t3 = reform(coeffmatrix[spindex.n4s,spindex.o2  , *])*zmajnow.o2den ; #
    n4s_loss_t4 = reform(coeffmatrix[spindex.n4s,spindex.no  , *])*zminornow.no ; #
    n4s_loss_t5 = reform(coeffmatrix[spindex.n4s,spindex.o2d_p, *])*zionnow.o2d_p
    n4s_loss_t6 = zpid.i_n4s
    
    n4s_prod = n4s_prod_t1  + n4s_prod_t2  + n4s_prod_t3  + $
               n4s_prod_t4  + n4s_prod_t5  + n4s_prod_t6  + $
               n4s_prod_t7  + n4s_prod_t8  + n4s_prod_t9  + $
               n4s_prod_t10 + n4s_prod_t11 + n4s_prod_t12 + $
               n4s_prod_t13 + n4s_prod_t14 + n4s_prod_t15 + $
               n4s_prod_t16 + n4s_prod_t17 + n4s_prod_t18 + $
               n4s_prod_t19 + n4s_prod_t20 + n4s_prod_t21

    n4s_loss = n4s_loss_t1 + n4s_loss_t2 + n4s_loss_t3 + n4s_loss_t4 + $
               n4s_loss_t5 + n4s_loss_t6 
    zminor.n4s_pce = ace_1d_pcecalc(zminornow.n4s_pce,n4s_prod,n4s_loss,model.del_time) ; <1e9
    
;    if(i EQ 3000) THEN STOP
    
    ;n_4s_temp = ace_1d_pcecalc(zminornow.n4s,n4s_prod,n4s_loss,model.del_time) ; needed for boundary condition
;   N4S is assumed to be in PCE equilibrium at lbc
;   n4s_q = 1.
;   Specify either lb1-lb4 or a flux
    ; n4s_lbc = zminor.n4s_pce[0] ; rhs
    ; n4s_lbc_flux = 0.;1e-7 ; If flux is non zero, it WILL be used in place of any other LBC
    ;zminor.n4s = ace_1d_minor_solv(zminornow.n4s,mass.n4s,zmajnow,n4s_prod,n4s_loss,k_m.n4s,flux.n4s, model, mass, pconst,i)
    ;zminor.n4s_mmr = ace_1d_cm3_mmr(zminor.n4s, mass.n4s, zmajnow, model, pconst)
    zminor.n4s_mmr = ace_1d_minor_n4s(zminornow, zmajnow, n4s_prod, n4s_loss, model, mass, pconst)
    zminor.n4s = ace_1d_mmr_cm3(zminor.n4s_mmr, mass.n4s, zmajnow, model, pconst)
    ;zminor.n4s = zminornow.n4s ; zminor.n4s_pce ; 

;;;;; N(2D)

    n2d_prod_t1 = 2.*zei.ped_n2*zmajnow.n2den*yieldmatrix[spindex.n2, spindex.pe, spindex.n2d]
    n2d_prod_t2 = 2.*zei.aed_n2*zmajnow.n2den*yieldmatrix[spindex.n2, spindex.ae, spindex.n2d]; check this yield 
    n2d_prod_t3 = 2.*zmajnow.n2den*zpid.j_n2*yieldmatrix[spindex.n2, spindex.ph, spindex.n2d]
    n2d_prod_t4 = zmajnow.n2den*zpid.di_n2*0.6
    n2d_prod_t5 = zmajnow.n2den*zei.di_n2*0.6
    
    n2d_prod_t6 = reform(ratematrix[spindex.n2_p,spindex.o, *])*yieldmatrix[spindex.n2_p,spindex.o,spindex.n2d] ;# 
    n2d_prod_t7 = reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p, spindex.o2, spindex.n2d] ; #
    n2d_prod_t8 = reform(ratematrix[spindex.n2p,spindex.o, *])*yieldmatrix[spindex.n2p, spindex.o, spindex.n2d] ; #      
    n2d_prod_t9 = reform(ratematrix[spindex.n2_p,spindex.e, *])*yieldmatrix[spindex.n2_p,spindex.e,spindex.n2d] ;#
    n2d_prod_t10 = reform(ratematrix[spindex.n2p,spindex.onebody, *])*yieldmatrix[spindex.n2p, spindex.onebody, spindex.n2d] ; #
    n2d_prod_t11 = reform(ratematrix[spindex.no_p,spindex.e, *])*yieldmatrix[spindex.no_p, spindex.e, spindex.n2d] ; #
    n2d_prod_t12 = reform(ratematrix[spindex.n2p,spindex.e, *])*n2p_n2d_yield ; #  

    n2d_prod_t13_v3 = k_n2a_o[3]*zminornow.n2a[3,*]*zmajnow.oden*n2a_no_n2d_yield ; production of N(2D) from N2A(v=3) by reaction with O
    n2d_prod_t13_v4 = k_n2a_o[4]*zminornow.n2a[4,*]*zmajnow.oden*n2a_no_n2d_yield
    n2d_prod_t13_v5 = k_n2a_o[5]*zminornow.n2a[5,*]*zmajnow.oden*n2a_no_n2d_yield
    n2d_prod_t13_v6 = k_n2a_o[6]*zminornow.n2a[6,*]*zmajnow.oden*n2a_no_n2d_yield

    
    n2d_loss_t1 = reform(coeffmatrix[spindex.n2d,spindex.o2, *])*zmajnow.o2den ; #
    n2d_loss_t2 = reform(coeffmatrix[spindex.n2d,spindex.e, *])*zionnow.e ; #      
    n2d_loss_t3 = reform(coeffmatrix[spindex.n2d,spindex.o, *])*zmajnow.oden ; #      
    n2d_loss_t4 = reform(coeffmatrix[spindex.n2d,spindex.no, *])*zminornow.no ; #
    n2d_loss_t5 = reform(coeffmatrix[spindex.n2d,spindex.onebody, *]) ; #  
    n2d_loss_t6 = reform(coeffmatrix[spindex.n2d,spindex.n2, *])*zmajnow.n2den ; # 
    n2d_loss_t7 = reform(coeffmatrix[spindex.n2d,spindex.o_p, *])*zionnow.o_p
    n2d_loss_t8 = reform(coeffmatrix[spindex.n2d,spindex.o2_p, *])*zionnow.o2_p

    n2d_prod = n2d_prod_t1 + n2d_prod_t2  + n2d_prod_t3  + n2d_prod_t4 + $
               n2d_prod_t5 + n2d_prod_t6  + n2d_prod_t7  + n2d_prod_t8 + $
               n2d_prod_t9 + n2d_prod_t10 + n2d_prod_t11 + n2d_prod_t12 + $
               n2d_prod_t13_v3 + n2d_prod_t13_v4 + n2d_prod_t13_v5 + n2d_prod_t13_v6

    n2d_loss = n2d_loss_t1 + n2d_loss_t2 + n2d_loss_t3 + $
               n2d_loss_t4 + n2d_loss_t5 + n2d_loss_t6 + $
               n2d_loss_t7 + n2d_loss_t8
 
 
    zminor.n2d = ace_1d_pcecalc(zminornow.n2d,n2d_prod,n2d_loss,model.del_time) ; fltarr(model.nlev) ;
;;;;; N(2P)

    n2p_prod_t1 = 2.*zmajnow.n2den*zei.ped_n2*yieldmatrix[spindex.n2, spindex.pe, spindex.n2p]
    n2p_prod_t2 = 2.*zmajnow.n2den*zpid.j_n2*yieldmatrix[spindex.n2, spindex.ph, spindex.n2p]
    n2p_prod_t3 = 2.*zmajnow.n2den*zei.aed_n2*yieldmatrix[spindex.n2, spindex.ae, spindex.n2p]
    
    ; n2p_prod_t4 = .25*zmajnow.n2den*zpid.di_n2 ; Need to check this yield before including
    ; n2p_prod_t5 = zei.di_n2*zmajnow.n2den ; Need to check this yield before including
    ; The N(2P) yield for this recombination is zero
    n2p_prod_t6 = reform(ratematrix[spindex.n2_p, spindex.e, *])*yieldmatrix[spindex.n2_p, spindex.e, spindex.n2p];#  

    n2p_loss_t1 = reform(coeffmatrix[spindex.n2p,spindex.o, *])*zmajnow.oden ; #
    n2p_loss_t2 = reform(coeffmatrix[spindex.n2p,spindex.o2, *])*zmajnow.o2den ; #
    n2p_loss_t3 = reform(coeffmatrix[spindex.n2p,spindex.n2, *])*zmajnow.n2den ; # FIND REFERENCE FOR THIS ONE   
    n2p_loss_t4 = reform(coeffmatrix[spindex.n2p,spindex.no, *])*zminornow.no  ; # FIND REFERENCE FOR THIS ONE
    n2p_loss_t5 = reform(coeffmatrix[spindex.n2p,spindex.e, *])*zionnow.e ; #      
    n2p_loss_t6 = reform(coeffmatrix[spindex.n2p,spindex.onebody, *]) ; #
    n2p_loss_t7 = reform(coeffmatrix[spindex.n2p,spindex.o2_p, *])*zionnow.o2_p ; #    
    
    n2p_prod = n2p_prod_t1 + n2p_prod_t2 + n2p_prod_t3 + n2p_prod_t6 ;+ n2p_prod_t4 + $
               ; n2p_prod_t5 

    n2p_loss = n2p_loss_t1 + n2p_loss_t2 + n2p_loss_t3 + n2p_loss_t4 + $
               n2p_loss_t5 + n2p_loss_t6 + n2p_loss_t7 

    zminor.n2p = ace_1d_pcecalc(zminornow.n2p,n2p_prod,n2p_loss,model.del_time)

END
