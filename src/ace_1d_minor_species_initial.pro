; Establish densities of minor species assuming photochemical equilibrium

; Minor species tracked in the model are 
; Nitric Oxide (NO)
; Atomic nitrogen
; - Ground State (N4S)
; - Electronically Excited (N(2D), N(2P))
; Electronically excited molecular nitrogen (N2(A))
; Electronically excited atomic nitrogen (O(1D))
; Carbon Dioxide (CO2)

; Nitric Oxide (NO)
k_n2a_o = [2.8, 3.4, 3.34, 3.48, 2.08, 4.26, 4.58, 7.4]*1e-11 ; Level specific reaction rate for v = 3-6
n2a_no_n2d_yield = 0.57
	
no_prod_t1 = reform(ratematrix[spindex.n2p,spindex.o2, *])
no_prod_t2 = reform(ratematrix[spindex.n2d,spindex.o2, *])
no_prod_t3 = reform(ratematrix[spindex.n4s,spindex.o2, *])
no_prod_t4_v3 = k_n2a_o[3]*zminor.n2a[3,*]*zmaj.oden*n2a_no_n2d_yield
no_prod_t4_v4 = k_n2a_o[4]*zminor.n2a[4,*]*zmaj.oden*n2a_no_n2d_yield
no_prod_t4_v5 = k_n2a_o[5]*zminor.n2a[5,*]*zmaj.oden*n2a_no_n2d_yield
no_prod_t4_v6 = k_n2a_o[6]*zminor.n2a[6,*]*zmaj.oden*n2a_no_n2d_yield
no_prod_t5 = reform(ratematrix[spindex.n_p,spindex.o2, *])*yieldmatrix[spindex.n_p,spindex.o2,spindex.no]

no_loss_t1 = zpid.i_no
no_loss_t2 = zpid.j_no
no_loss_t3 = reform(coeffmatrix[spindex.no,spindex.n2_p, *])*zion.n2_p
no_loss_t4 = reform(coeffmatrix[spindex.no,spindex.n_p , *])*zion.n_p 
no_loss_t5 = reform(coeffmatrix[spindex.no,spindex.o_p , *])*zion.o_p 
no_loss_t6 = reform(coeffmatrix[spindex.no,spindex.o2_p, *])*zion.o2_p  
no_loss_t7 = reform(coeffmatrix[spindex.no,spindex.n2p , *])*zminor.n2p 
no_loss_t8 = reform(coeffmatrix[spindex.no,spindex.n2d , *])*zminor.n2d 
no_loss_t9 = reform(coeffmatrix[spindex.no,spindex.n4s , *])*zminor.n4s 
no_loss_t10 = reform(coeffmatrix[spindex.no,spindex.o2d_p , *])*zion.o2d_p

no_prod = no_prod_t1 + no_prod_t2 + no_prod_t3 + no_prod_t4_v3 $
+ no_prod_t4_v4 + no_prod_t4_v5 + no_prod_t4_v6 + no_prod_t5

no_loss = no_loss_t1 + no_loss_t2 + no_loss_t3 + no_loss_t4 + $
no_loss_t5 + no_loss_t6 + no_loss_t7 + no_loss_t8 + $
no_loss_t9 + no_loss_t10

zminor.no_pce = ace_1d_pcecalc(zminor.no_pce, no_prod, no_loss, model.del_time)


;;; N(4S)

    n2p_n4s_yield = 1.6e-12*(zion.te)^0.85/(1.6e-12*(zion.te)^0.85 + 9.5e-9) > 0.
    n2p_n2d_yield = 9.5e-9/(1.6e-12*(zion.te)^0.85 + 9.5e-9) > 0.
    
    n4s_prod_t1 = 2.*zei.ped_n2*zmaj.n2den*yieldmatrix[spindex.n2, spindex.pe, spindex.n4s]
    n4s_prod_t2 = 2.*zei.aed_n2*zmaj.n2den*yieldmatrix[spindex.n2, spindex.ae, spindex.n4s]
    n4s_prod_t3 = zmaj.n2den*zpid.di_n2*0.4 ; Photoelectron dissociative ionization of N2, times N(4S) branching ratio
    n4s_prod_t4 = zmaj.n2den*zei.di_n2*0.4 ; Dissociative ionization of N2 by EUV photons, times N(4S) branching ratio
    n4s_prod_t5 = 2.*zmaj.n2den*zpid.j_n2*yieldmatrix[spindex.n2, spindex.ph, spindex.n4s] ; Photodissociation of N2
    
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
    n4s_prod_t20 = zminor.no*zpid.j_no ; NO photolysis    
    n4s_prod_t21 = reform(ratematrix[spindex.n2p,spindex.o2_p, *]) ; #
    
    n4s_loss_t1 = reform(coeffmatrix[spindex.n4s,spindex.n2_p, *])*zion.n2_p ;#
    n4s_loss_t2 = reform(coeffmatrix[spindex.n4s,spindex.o2_p, *])*zion.o2_p ; #    
    n4s_loss_t3 = reform(coeffmatrix[spindex.n4s,spindex.o2  , *])*zmaj.o2den ; #
    n4s_loss_t4 = reform(coeffmatrix[spindex.n4s,spindex.no  , *])*zminor.no ; #
    n4s_loss_t5 = reform(coeffmatrix[spindex.n4s,spindex.o2d_p, *])*zion.o2d_p
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
    zminor.n4s_pce = ace_1d_pcecalc(zminor.n4s_pce,n4s_prod,n4s_loss,model.del_time) ; <1e9




STOP
