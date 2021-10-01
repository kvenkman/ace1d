PRO ace_1d_updatechem, zmaj, zion, zminor, model,  $
                       spindex, exomatrix, coeffmatrix, ratematrix, $
                       heatmatrix
                       
; Reaction coefficients, adopted from 1D NO model as of March, 2016
; Scott, Karthik, and Justin Y

; Conventions for adding reactions to the model
; ------------------------------------------
; Each reaction block begins writing out the reaction
; Check that the reaction isn't already included in the model
; Then, reference the source from which the rate coefficient is obtained
; Any other comments are added after this
; Define the rate coefficient in the variable "coeff"
; Pass coeff along with the reactant and product indices to ace_1d_updatecrh

coeffmatrix = fltarr(spindex.nsp,spindex.nsp,model.nlev)
ratematrix = fltarr(spindex.nsp,spindex.nsp,model.nlev)
heatmatrix = fltarr(spindex.nsp,spindex.nsp,model.nlev)

;  Ion Neutral Reactions
;  ---------------------
;  
; N+ + NO -> N(4S) + NO+ / N2+ + O / N2 + O+
; Yonker [2013], Table D4
coeff = 6.5e-9*(zion.ti)^(-.44)
ace_1d_update_crh, spindex.n_p, spindex.no, zion.n_p, zminor.no, $
                   coeff, exomatrix[spindex.n_p, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; N+ + O(3P) -> O+ + N(4S)
; Yonker [2013], Table D4
coeff = 4.5e-12
ace_1d_update_crh, spindex.n_p, spindex.o, zion.n_p, zmaj.oden, $
                   coeff, exomatrix[spindex.n_p, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix

; N+ + O2 -> N(2D) + O2+ / N(4S) + O2+ / NO+ + O / NO + O+(4S)
; Yonker [2013]
coeff = 5.5e-10
ace_1d_update_crh, spindex.n_p, spindex.o2, zion.n_p, zmaj.o2den, $
                   coeff, exomatrix[spindex.n_p, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

; N2+ + e -> N(2D) + N(2D) / N(2D) + N(4S) / N(2P) + N(4S)
; Yonker [2013]
coeff = 2.2e-7*(300./zion.te)^(0.39)
ace_1d_update_crh, spindex.n2_p, spindex.e, zion.n2_p, zion.e, $
                   coeff, exomatrix[spindex.n2_p, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; N2+ + N(4S) -> N+ + N2
; Yonker [2013]
coeff = 1e-11
ace_1d_update_crh, spindex.n2_p, spindex.n4s, zion.n2_p, zminor.n4s, $
                   coeff, exomatrix[spindex.n2_p, spindex.n4s], $
                   coeffmatrix, ratematrix, heatmatrix

; N2+ + NO -> NO+ + N2       
; Yonker [2013]
coeff = 7.5e-9*(zion.ti)^(-.52)
ace_1d_update_crh, spindex.n2_p, spindex.no, zion.n2_p, zminor.no, $
                   coeff, exomatrix[spindex.n2_p, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; N2+ + O(3P) -> NO+ + N(2D) ;
; Yonker [2013], expanded using Fox & Sung [2001]
coeff = fltarr(model.nlev)
mask1 = (zion.ti LE 1500)
mask2 = (zion.ti GT 1500)
coeff = mask1*1.33e-10*(300./zion.ti)^(0.44) + $
        mask2*6.55e-11*(zion.ti/300.)^(0.2)
ace_1d_update_crh, spindex.n2_p, spindex.o, zion.n2_p, zmaj.oden, $
                   coeff, exomatrix[spindex.n2_p, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix

; N2+ + O2 -> O2+ + N2 ;
; Expanded from Yonker [2013] using P. G. Richards [2011]
coeff = fltarr(model.nlev)
mask1 = zion.ti LE 1000
mask2 = zion.ti GT 1000
coeff = (mask1*5.1e-11*(300./zion.ti)^(1.16) + $
        mask2*1.26e-11*(zion.ti/1000.)^0.67)
;coeff = 5.1e-11*(300./zion.ti)^(1.16)
ace_1d_update_crh, spindex.n2_p, spindex.o2, zion.n2_p, zmaj.o2den, $
                   coeff, exomatrix[spindex.n2_p, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

; NO+ + e -> O + N(2D)/N(4S)
; Yonker [2013]
coeff = 3.5e-7*(300./zion.te)^(0.69)
ace_1d_update_crh, spindex.no_p, spindex.e, zion.no_p, zion.e, $
                   coeff, exomatrix[spindex.no_p, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; O+ + N(2D) -> N+ + O
; P. G. Richards [2011], Table 1
coeff = 1.3e-10 ; fltarr(model.nlev); 
ace_1d_update_crh, spindex.o_p, spindex.n2d, zion.o_p, zminor.n2d, $
                   coeff, exomatrix[spindex.o_p, spindex.n2d], $
                   coeffmatrix, ratematrix, heatmatrix

; O+ + N2 -> NO+ + N(4S)                     
; Expanded from Yonker [2013] using P. G. Richards [2011]
mask1 = zion.ti LE 1000
mask2 = zion.ti GT 1000
coeff = fltarr(model.nlev)
coeff = (mask1*1.2e-12*(300./zion.ti)^0.45 + $
        mask2*7e-13*(zion.ti/1000.)^2.12)
;coeff = 1.2e-12*(300./zion.ti)^0.45
ace_1d_update_crh, spindex.o_p, spindex.n2, zion.o_p, zmaj.n2den, $
                   coeff, exomatrix[spindex.o_p, spindex.n2], $
                   coeffmatrix, ratematrix, heatmatrix

; O+ + NO -> NO+ + O
; Yonker [2013]
;coeff = fltarr(model.nlev);
coeff = 5.01e-13*(300./zion.ti)^1.68; + 4.02e-12*exp(-901./zion.ti)
ace_1d_update_crh, spindex.o_p, spindex.no, zion.o_p, zminor.no, $
                   coeff, exomatrix[spindex.o_p, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; O+ + O2 -> O2+ + O
; Yonker [2013]
coeff = fltarr(model.nlev)
coeff = 1.7e-11*(300./zion.ti)^0.77 + 8.54e-11*exp(-3461./zion.ti)
ace_1d_update_crh, spindex.o_p, spindex.o2, zion.o_p, zmaj.o2den, $
                   coeff, exomatrix[spindex.o_p, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2D) + e -> O+(4S) + e ;
; P. G. Richards [2011], Table 1
coeff  = 6.03e-8*sqrt(300./zion.te)
ace_1d_update_crh, spindex.o2d_p, spindex.e, zion.o2d_p, zion.e, $
                   coeff, exomatrix[spindex.o2d_p, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2D) + N(4S) -> N+ + O
; P. G. Richards [2011], Table 1
coeff  = 1.5e-10
ace_1d_update_crh, spindex.o2d_p, spindex.n4s, zion.o2d_p, zminor.n4s, $
                   coeff, exomatrix[spindex.o2d_p, spindex.n4s], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2D) + N2 -> O(3P) + N2+ 
; see Fox and Sung (2001)
coeff = 5.7e-10*exp(-400./zion.ti)
ace_1d_update_crh, spindex.o2d_p, spindex.n2, zion.o2d_p, zmaj.n2den, $
                   coeff, exomatrix[spindex.o2d_p, spindex.n2], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2D) + NO -> NO+ + O
; P. G. Richards [2011], Table 1
coeff  = 1.2e-9 
ace_1d_update_crh, spindex.o2d_p, spindex.no, zion.o2d_p, zminor.no, $
                   coeff, exomatrix[spindex.o2d_p, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2D) + O(3P) -> O+(4S) + O(3P) ;
; P. G. Richards [2011] ; Yonker [2013]
coeff  = 1.0e-11
ace_1d_update_crh, spindex.o2d_p, spindex.o, zion.o2d_p, zmaj.oden, $
                   coeff, exomatrix[spindex.o2d_p, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2D) + O2 -> O2+ + O  
; from Johnsen and Biondi, JCP, (1980)
coeff = 7.e-10
ace_1d_update_crh, spindex.o2d_p, spindex.o2, zion.o2d_p, zmaj.o2den, $
                   coeff, exomatrix[spindex.o2d_p, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2P) + e -> O+(4S) + e  / O+(2D) + e
; P. G. Richards [2011], Table 1
coeff  = (1.84e-7*(300./zion.te)^0.5 + 3.03e-8*(300./zion.te)^0.5)
ace_1d_update_crh, spindex.o2p_p, spindex.e, zion.o2p_p, zion.e, $
                   coeff, exomatrix[spindex.o2p_p, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2P) + N2 -> O(3P) + N2+ 
; see Fox and Sung (2001)
coeff = 5.7e-10*exp(-400./zion.ti)
ace_1d_update_crh, spindex.o2p_p, spindex.n2, zion.o2p_p, zmaj.n2den, $
                   coeff, exomatrix[spindex.o2p_p, spindex.n2], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2P) + O(3P) -> O+(4S) + O(3P) ;
; Yonker [2013]
coeff  = 5.0e-11
ace_1d_update_crh, spindex.o2p_p, spindex.o, zion.o2p_p, zmaj.oden, $
                   coeff, exomatrix[spindex.o2p_p, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2P) + O2 -> O2+ + O  
; from Johnsen and Biondi, JCP, (1980)
coeff = 7.e-10
ace_1d_update_crh, spindex.o2p_p, spindex.o2, zion.o2p_p, zmaj.o2den, $
                   coeff, exomatrix[spindex.o2p_p, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

; O2+ + e -> O + O
; Fennelly and Torr [1994]
coeff = fltarr(model.nlev)
mask1 = zion.te LT 1200
mask2 = zion.te GE 1200
coeff = mask1*(2e-7*(300./zion.te)^0.7) + $
        mask2*(1.6e-7*(300./zion.te)^0.55)
ace_1d_update_crh, spindex.o2_p, spindex.e, zion.o2_p, zion.e, $
                   coeff, exomatrix[spindex.o2_p, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; O2+ + N(2D) -> NO+ + O / N+ + O2
; P. G. Richards [2011], Table 1 / Goldan [1966] / Fox and Sung [2001]
coeff = 1.8e-10 + 8.65e-11
ace_1d_update_crh, spindex.o2_p, spindex.n2d, zion.o2_p, zminor.n2d, $
                   coeff, exomatrix[spindex.o2_p, spindex.n2d], $
                   coeffmatrix, ratematrix, heatmatrix

; O2+ + N(2P) -> N + O2+
; P. G. Richards [2011], Table 1, Zipf et al [1980]
coeff = 2.2e-11
ace_1d_update_crh, spindex.o2_p, spindex.n2p, zion.o2_p, zminor.n2p, $
                   coeff, exomatrix[spindex.o2_p, spindex.n2p], $
                   coeffmatrix, ratematrix, heatmatrix

; O2+ + N(4S) -> NO+ + O
; Yonker [2013]
coeff = 1.33e-10
ace_1d_update_crh, spindex.o2_p, spindex.n4s, zion.o2_p, zminor.n4s, $
                   coeff, exomatrix[spindex.o2_p, spindex.n4s], $
                   coeffmatrix, ratematrix, heatmatrix

; O2+ + NO -> NO+ + O2
; Yonker [2013]
coeff = 4.5E-10
ace_1d_update_crh, spindex.o2_p, spindex.no, zion.o2_p, zminor.no, $
                   coeff, exomatrix[spindex.o2_p, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix
                   
; N2+ + e -> N2(A), 8 vibrational levels?
;chema.a_8 = [0.422106,0.377770,0.347252,0.323350,0.311022,0.305038,0.299031,0.422961]    
;exoa.a_8  = 0.0000000000000000000000000000
;rate = ;????????????????????????????????
;ace_1d_reaction,spindex.n2_p,spindex.e,spindex.n2a,spindex.e,rate,exoa.a_8,chemmat,heatmat,spname,reaction


; Neutral - Neutral Reactions
; ---------------------------
;
; N(2P) + O2 -> NO + O (3P,1D)
; Yonker [2013]
coeff = 3.1e-12*exp(-60./zmaj.tn)
ace_1d_update_crh, spindex.n2p, spindex.o2, zminor.n2p, zmaj.o2den, $
                   coeff, exomatrix[spindex.n2p, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

;N(2P) + O(3P) -> N(4S,2D) + O / NO+ + e
; Yonker [2013]
coeff = 2.7e-11
ace_1d_update_crh, spindex.n2p, spindex.o, zminor.n2p, zmaj.oden, $
                   coeff, exomatrix[spindex.n2p, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix

; N(2P) + e -> N(2D)/N(4S) + e
; Yonker [2013]  
coeff = 9.5e-9 + 1.6e-12*(zion.te)^0.85
ace_1d_update_crh, spindex.n2p, spindex.e, zminor.n2p, zion.e, $
                   coeff, exomatrix[spindex.n2p, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; N(2D) + N2 -> N(4S) + N2;
; Yonker [2013]
coeff = 1.74e-14
ace_1d_update_crh, spindex.n2d, spindex.n2, zminor.n2d, zmaj.n2den, $
                   coeff, exomatrix[spindex.n2d, spindex.n2], $
                   coeffmatrix, ratematrix, heatmatrix

; N(2D) + O2 -> NO + O(3P)
; Duff (2003) (t3)
coeff = 6.2e-12*(zmaj.tn/300.) ;  
;coeff = 5e-12
ace_1d_update_crh, spindex.n2d, spindex.o2, zminor.n2d, zmaj.o2den, $
                   coeff, exomatrix[spindex.n2d, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix


; N(2D) + O -> N(4S) + O   
; pg77 / Table D.6, JY dissertation the Herron recommendation for the
; temperature dependence is used, but fit to the Fell et al [1990]
;  value at room temperature
coeff = 1.65E-12*exp(-260./zmaj.tn)
; Roble 1995
;coeff = 4.5e-13
ace_1d_update_crh, spindex.n2d, spindex.o, zminor.n2d, zmaj.oden, $
                   coeff, exomatrix[spindex.n2d, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix

; N(2D) + NO -> N2 + O(3P)
; Yonker [2013]
coeff = 6.7E-11;  (t3)
ace_1d_update_crh, spindex.n2d, spindex.no, zminor.n2d, zminor.no, $
                   coeff, exomatrix[spindex.n2d, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; N(2D) + e -> N(4S) + e    
; Yonker [2013]
coeff = 3.86e-10*(zion.te/300.)^.81
; Roble 1995
;coeff = 3.6e-10*(zion.te/300.)^0.5
;coeff = 3.8e-12*(zion.te)^.81
ace_1d_update_crh, spindex.n2d, spindex.e, zminor.n2d, zion.e, $
                   coeff, exomatrix[spindex.n2d, spindex.e], $
                   coeffmatrix, ratematrix, heatmatrix

; N(4S) + O2 -> NO + O
; Sander et al. [2006]
coeff = 1.5E-11*exp(-3600./zmaj.tn)
ace_1d_update_crh, spindex.n4s, spindex.o2, zminor.n4s, zmaj.o2den, $
                   coeff, exomatrix[spindex.n4s, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix
       
; N(4S) + NO -> N2 + O(3P)
; this is the jpl 2006 recommendation
coeff = 2.1e-11*exp(100./zmaj.tn)
ace_1d_update_crh, spindex.n4s, spindex.no, zminor.n4s, zminor.no, $
                   coeff, exomatrix[spindex.n4s, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; O + O + M -> O2 + M
; M is N2 or O2. Roble [1987]
; This is defined separately as it is a 3 body reaction
coeff = 9.59e-34*exp(480./zmaj.tn)
coeffmatrix[spindex.o, spindex.o, *] = coeff
ratematrix[spindex.o, spindex.o, *] = coeff*zmaj.oden*zmaj.oden*(zmaj.n2den+zmaj.o2den)
heatmatrix[spindex.o,spindex.o, *] = exomatrix[spindex.o,spindex.o]*ratematrix[spindex.o, spindex.o, *]


; This reaction isn't important for the current version of the model
; O + O2 + M -> O3 + M
; Herron, with t-dep
; coeff = 6.0e-34 * (zmaj.tn/300.)^(-2.4)
; ace_1d_update_crh, spindex.o, spindex.o2, zmaj.oden, zmaj.o2den, $
;                  coeff, exomatrix[spindex.o, spindex.o2], $
;                  coeffmatrix, ratematrix, heatmatrix

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; N2(A) reactions ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; N2(A) loss to O(3P), O2 and radiative losses
; N2(A) + O -> NO + N(2D)
; The rate coeﬃcents are shown in Tables D.6, C.1, and D.5. of Yonker's Dissertation
 coeff = 1e-13
 ace_1d_update_crh, spindex.n2a, spindex.o, zminor.n2a, zmaj.oden, $
                   coeff, exomatrix[spindex.n2a, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix


; Quenching Reactions
; -------------------

; note that quenching reactions have a reactant and product that is the same - 
; need to make sure this is handled right in P and L calculations


; N(2P) + N2    will never be important 
; Herron 1999
coeff = 5e-17
ace_1d_update_crh, spindex.n2p, spindex.n2, zminor.n2p, zmaj.n2den, $
                   coeff, exomatrix[spindex.n2p, spindex.n2], $
                   coeffmatrix, ratematrix, heatmatrix

; N(2P) + NO
; Herron 1999
coeff = 2.9e-11
ace_1d_update_crh, spindex.n2p, spindex.no, zminor.n2p, zminor.no, $
                   coeff, exomatrix[spindex.n2p, spindex.no], $
                   coeffmatrix, ratematrix, heatmatrix

; O(1D) quenching rate coefficients are obtained from Fennelly and Torr [1994] / Roble [1995]
; O(1D) + N2
coeff =  2e-11*exp(107.8/zmaj.tn)
ace_1d_update_crh, spindex.o1d, spindex.n2, zminor.o1d, zmaj.n2den, $
                   coeff, exomatrix[spindex.o1d, spindex.n2], $
                   coeffmatrix, ratematrix, heatmatrix
                   

; O(1D) + O2
coeff =  2.9e-11*exp(67.5/zmaj.tn)
ace_1d_update_crh, spindex.o1d, spindex.o2, zminor.o1d, zmaj.o2den, $
                   coeff, exomatrix[spindex.o1d, spindex.o2], $
                   coeffmatrix, ratematrix, heatmatrix

; O(1D) + O
coeff =  8e-12
ace_1d_update_crh, spindex.o1d, spindex.o, zminor.o1d, zmaj.oden, $
                   coeff, exomatrix[spindex.o1d, spindex.o], $
                   coeffmatrix, ratematrix, heatmatrix
                   
; Radiative Relaxation
; ------------------------

; N(2D) -> N(4S) + hv   
coeff = (fltarr(model.nlev) + 1.279e-5)
ace_1d_update_crh, spindex.n2d, spindex.onebody, zminor.n2d, zminor.n2d, $
                   coeff, exomatrix[spindex.n2d, spindex.onebody], $
                   coeffmatrix, ratematrix, heatmatrix

; The yield matrix will take care of the N(2D)/N(4S) channels
; N(2P) -> N(2D)/N(4S) + hv
coeff = (fltarr(model.nlev) + 8.054e-2 + 5.308e-3)
ace_1d_update_crh, spindex.n2p, spindex.onebody, zminor.n2p, zminor.n2p, $
                   coeff, exomatrix[spindex.n2p, spindex.onebody], $
                   coeffmatrix, ratematrix, heatmatrix

; O+(2P) -> O+(4S)/O+(2D) + hv
; Roble [1995]
coeff = (fltarr(model.nlev) + 0.047 + 0.171)
ace_1d_update_crh, spindex.o2p_p, spindex.onebody, zion.o2p_p, zion.o2p_p, $
                   coeff, exomatrix[spindex.o2p_p, spindex.onebody], $
                   coeffmatrix, ratematrix, heatmatrix
                   
; O(1D) -> O(3P) + hv (630 / 636.4 / 639.2)
; Fennelly and Torr [1994]
coeff = 0.0059 + 0.0018
ace_1d_update_crh, spindex.o1d, spindex.onebody, zminor.o1d, zminor.o1d, $
                   coeff, exomatrix[spindex.o1d, spindex.onebody], $
                   coeffmatrix, ratematrix, heatmatrix

END

; The following is mesosphere chemistry. This needs to be made to follow the same
; naming convenstions as above.
;;;;;New addtions 6/2013 Ox, HOx, and NOx for the mesosphere....
; reaction rates are taken from Brasseur and Solomon, 2005
;g_11_1_M = 6.0e-34 *((zmaj.tn/300.)^(-2.4)) 
;g_11_12 = 8.0e-12 * exp(-2060./zmaj.tn)
;g_o2_o2 = 2.22e-18 * ((maj.tn/300.)^0.78 ) ; this and the above from Thomas et al. 1984
;g_10_2 =1.8e-11* exp(110./zmaj.tn); quenching of o1d by n2, produces O(3P) + N2
;g_10_1 = 3.2e-11* exp(70./zmaj.tn) ; quenching of o1d by o2, produces O(3P) + O21Sigma
;g_10_12 = 1.2e-10 ; O(1d) + O3 -> O2 + O2

;g_11_11_m = 9.59e-34*exp(480./zmaj.tn)
;g_o21d_o23p =  0.0 ; deactivation by collision
;g_o1d_o3p = 0.0; quenching of o1d by o, produces O2
;a1_27 = 2.58e-4; radiative relaxation of O21Delta to the O23Sigma (ground state)

;g_h2o_o1d=2.2e-10 ; H2O + O(1D) -> 2OH
;g_h_o2_m=5.7e-32*((zmaj.tn/300.)^(-1.6)) ; H+O2+M -> HO2
;g_h_o3=1.4e-10*exp(-470./zmaj.tn) ; H + O3 -> OH + O2
;g_oh_o=2.2e-11*exp(120./zmaj.tn) ; OH + O -> O2 + H
;g_oh_o3=1.7e-12*exp(-940./zmaj.tn) ; OH + O3 -> O2 + HO2
;g_ho2_o3=1.0e-14*exp(-490./zmaj.tn) ; HO2 + O3 -> 2O2 + OH
;g_ho2_no=3.5e-12*exp(250./zmaj.tn) ; HO2 + NO -> NO2 + OH
;g_h_ho2=8.1e-11 ; H + HO2  -> ...
;g_oh_ho2=4.8e-11*exp(250./zmaj.tn) ; OH + HO2 -> H2O + O2
;g_oh_oh=4.2e-12*exp(-240./zmaj.tn) ; OH + OH -> H2O + O
;g_ho2_o=3.0e-11*exp(200./zmaj.tn)

;g_no_o3=2e-12*exp(-1400/zmaj.tn)
;g_no2_o=6.5e-12*exp(120./zmaj.tn)
;g_no2_o3=1.2e-13*exp(-2450./zmaj.tn)
;g_n_no2=5.8e-12*exp(220./zmaj.tn)
;g_n4s_oh=5.5e-11

;        g_11_1_m  : g_11_1_m ,$
;        g_11_12   : g_11_12  ,$
;        g_10_2    : g_10_2   ,$
;        g_10_1    : g_10_1   ,$
;        g_10_12   : g_10_12  ,$
;        g_11_11_m : g_11_11_m,$
;        g_h2o_o1d : g_h2o_o1d,$
;        g_h_o2_m  : g_h_o2_m ,$
;        g_h_o3    : g_h_o3   ,$
;        g_oh_o    : g_oh_o   ,$
;        g_oh_o3   : g_oh_o3  ,$
;        g_ho2_o3  : g_ho2_o3 ,$
;        g_ho2_no  : g_ho2_no ,$
;        g_h_ho2   : g_h_ho2  ,$
;        g_oh_ho2  : g_oh_ho2 ,$
;        g_oh_oh   : g_oh_oh  ,$
;        g_ho2_o   : g_ho2_o  ,$
;        g_no_o3   : g_no_o3  ,$
;        g_no2_o   : g_no2_o  ,$
;        g_no2_o3  : g_no2_o3 ,$
;        g_n_no2   : g_n_no2  ,$
;        g_n4s_oh  : g_n4s_oh  $
