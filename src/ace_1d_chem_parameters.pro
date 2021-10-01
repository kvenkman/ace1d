; Define yields and exothermicities of all reactions utilized in the model
; This code is run once, at model initialization

; Channel yields are dimensionless; exothermicities are in eV

;  Ion Neutral Reactions
; -------------------------
;;;;;;;;;; N+ reactions ;;;;;;;;;;;;;;;
; N+ + NO -> NO+ + N
; new reaction added 9/2008 by jy.  reference Midey, Miller, Viggiano, JCP, 2004
yield1 = 0.91 ; NO+ + N(4S) channel, 5.3eV
yield2 = 0.07 ; O + N2+ channel, 2.2eV
yield3 = 1. - yield1 - yield2 ; N2 + O+(4S) channel, 4.2eV

ace_1d_yieldmatrix, spindex.n_p, spindex.no, spindex.no_p, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.no, spindex.n4s, yield1, yieldmatrix

ace_1d_yieldmatrix, spindex.n_p, spindex.no, spindex.o, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.no, spindex.n2_p, yield2, yieldmatrix

ace_1d_yieldmatrix, spindex.n_p, spindex.no, spindex.n2, yield3, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.no, spindex.o_p, yield3, yieldmatrix

exo1 = 5.3
exo2 = 2.2
exo3 = 4.2

exo = yield1*exo1 + yield2*exo2 + yield3*exo3
ace_1d_exo_matrix, spindex.n_p, spindex.no, exo, exomatrix

; N+ + O -> O+ + N(4S)
; Roble, 1995
exo = 0.98
ace_1d_exo_matrix, spindex.n_p, spindex.o, exo, exomatrix

; N+ + O2 -> NO+ + O                  
; Yields have been taken from Yonker [2013]
yield1 = 0.5 ; into the N(2D) + O2+ channel.
yield2 = 0.42 ; into the O + NO+ channel
yield3 = 1. - yield1 - yield2 ; NO + O+(4S) channel
ace_1d_yieldmatrix, spindex.n_p, spindex.o2, spindex.n2d, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.o2, spindex.o2_p, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.o2, spindex.o, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.o2, spindex.no_p, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.o2, spindex.no, yield3, yieldmatrix
ace_1d_yieldmatrix, spindex.n_p, spindex.o2, spindex.o_p, yield3, yieldmatrix 

; From Midey and Viggiano [2006]
exo1 = 0.04 ; N(2D) + O2+ channel.
exo2 = 6.7 ; O + NO+ channel
exo3 = 2.3 ;  NO + O+(4S) channel

exo = yield1*exo1 + yield2*exo2 + yield3*exo3
ace_1d_exo_matrix, spindex.n_p, spindex.o2, exo, exomatrix

;;;;;;;;;; N2+ reactions ;;;;;;;;;;;;;;;

;N2+ + e -> N(2D) + N(4S)         
yield1 = 0.52 ; N(2D) + N(2D) channel
yield2 = 0.37 ; N(2D) + N(4S) channel
yield3 = 1. - yield1 - yield2 ; N(4S) + N(2P) channel
ace_1d_yieldmatrix, spindex.n2_p, spindex.e, spindex.n2d, 2.*yield1 + yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2_p, spindex.e, spindex.n4s, yield2 + yield3, yieldmatrix
ace_1d_yieldmatrix, spindex.n2_p, spindex.e, spindex.n2p, yield3, yieldmatrix

exo1 = 3.44
exo2 = 5.82
exo3 = 4.63
exo = exo1*yield1 + exo2*yield2 + exo3*yield3
ace_1d_exo_matrix, spindex.n2_p, spindex.e, exo, exomatrix

; N2+ + N -> N+ + N2
; Fox and Sung (2001)
; Value presently unknown
exo = 0.00000000000000000000000000000000000000
ace_1d_exo_matrix, spindex.n2_p, spindex.n4s, exo, exomatrix

; N2+ + NO -> NO+ + N2
exo = 6.3 ; Midey [2004]
ace_1d_exo_matrix, spindex.n2_p, spindex.no, exo, exomatrix

; N2+ + O -> NO+ + N(2D) / NO+ + N(4S) / O+ + N2
yield1 = 0.9 ;  NO+ + N(2D)
yield2 = 0.05 ; NO+ + N(4S) channel
yield3 = 1. - yield1 - yield2 ; O+ + N2 channel, Scott et al [1999]
ace_1d_yieldmatrix, spindex.n2_p, spindex.o, spindex.no_p, yield1+yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2_p, spindex.o, spindex.n2d, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n2_p, spindex.o, spindex.n4s, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2_p, spindex.o, spindex.o_p, yield3, yieldmatrix
ace_1d_yieldmatrix, spindex.n2_p, spindex.o, spindex.n2, yield3, yieldmatrix

; These need to be re-verified. Calculated from Scott et al. [1999], Table 1. 
; Assumed E[N(2D)] - E[N(4S)] = 2.38 eV
exo1 = 0.646
exo2 = 3.06
exo3 = 1.02

exo  = yield1*exo1 + yield2*exo2 + yield3*exo3
ace_1d_exo_matrix, spindex.n2_p, spindex.o, exo, exomatrix

; N2+ + O2 -> O2+ + N2 ;
; Dotan et al [1997]
exo = 3.5
ace_1d_exo_matrix, spindex.n2_p, spindex.o2, exo, exomatrix

;;;;;;;;;; NO+ reactions ;;;;;;;;;;;
; NO+ + e -> O + N(2D)/N(4S)
; Sheehan and ST.Maurice(2004)
; Hellberg et al.
yield1 = 0.95 ; N(2D) channel
yield2 = 1. - yield1 ; N(4S) channel
ace_1d_yieldmatrix, spindex.no_p, spindex.e, spindex.n2d, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.no_p, spindex.e, spindex.n4s, yield2, yieldmatrix

; Roble [1995]
exo1 = 0.38
exo2 = 2.75
exo = yield1*exo1 + yield2*exo2
ace_1d_exo_matrix, spindex.no_p, spindex.e, exo, exomatrix

;;;;;;;;;; O+ reactions ;;;;;;;;;;;

; O+ + N(2D) -> N+ + O      
; This value is used in TIE-GCM, but needs a reference. Bates [1989] doesn't mention this
; Roble [1995]
exo  = 1.45
ace_1d_exo_matrix, spindex.o_p, spindex.n2d, exo, exomatrix

; O+ + N2 -> NO+ + N(4S)                     
; This value is used in TIE-GCM, but needs a reference
; Roble [1995]
exo  = 1.0888
ace_1d_exo_matrix, spindex.o_p, spindex.n2, exo, exomatrix

; O+ + NO -> NO+ + O
; Dotan and Viggiano, JCP, 1999
; new reaction added 9/2008 by jy.  
exo = 4.3
ace_1d_exo_matrix, spindex.o_p, spindex.no, exo, exomatrix

; O+ + O2 -> O2+ + O
; This value is used in TIE-GCM, but needs a reference
; Roble [1995]
exo = 1.556
ace_1d_exo_matrix, spindex.o_p, spindex.o2, exo, exomatrix

; O+(2D) + e -> O+(4S) + e + 3.31 ev
; Roble 95, table 1, reaction 26
exo = 3.31
ace_1d_exo_matrix, spindex.o2d_p, spindex.e, exo, exomatrix

; O+(2D) + N(4S) -> N+ + O
exo = 0.
ace_1d_exo_matrix, spindex.o2d_p, spindex.n4s, exo, exomatrix

; O+(2D)+ N2 -> N2+ + O
; TIE-GCM
; Roble 95, Eq. 23
exo = 1.35
ace_1d_exo_matrix, spindex.o2d_p, spindex.n2, exo, exomatrix

; O+(2D) + NO -> NO+ + O
exo = 0.
ace_1d_exo_matrix, spindex.o2d_p, spindex.no, exo, exomatrix

; O+(2D) + O(3P) -> O+(4S) + O(3P) ;
; Stephan et al 2003
exo = 3.31
ace_1d_exo_matrix, spindex.o2d_p, spindex.o, exo, exomatrix

; O+(2D) + O2 -> O2+ + O  
; Roble 95
exo = 4.865
ace_1d_exo_matrix, spindex.o2d_p, spindex.o2, exo, exomatrix

; O+(2P) + e -> O+(4S) / O+(2D)
yield1 = 0.85860
yield2 = 1. - yield1
ace_1d_yieldmatrix, spindex.o2p_p, spindex.e, spindex.o2d_p, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.o2p_p, spindex.e, spindex.o_p, yield2, yieldmatrix

exo = 5.*yield1 + 1.69*yield2 ; Stephan et al. [2003]
ace_1d_exo_matrix, spindex.o2p_p, spindex.e, exo, exomatrix

; O+(2P) + O(3P) -> O+(4S) + O(3P) ;
; Stephan et al 2003
exo = 5.
ace_1d_exo_matrix, spindex.o2p_p, spindex.o, exo, exomatrix

; O+(2P) + O2 -> O2+ + O  
exo = 0.
ace_1d_exo_matrix, spindex.o2p_p, spindex.o2, exo, exomatrix

; O+(2P) + N2 -> N2+ + O
; Fox & Sung [2001] & reference therein (Li et. al [1997])
exo = 3.05
ace_1d_exo_matrix, spindex.o2p_p, spindex.n2, exo, exomatrix

;;;;;;;;;; O2+ reactions ;;;;;;;;;;;

; O2+ + e -> O(3P) + O(3P) / O(1D) + O(3P) / O(1D) + O(1D) 
; Schunk [2000], section 8.4. 
; We count the O(1S) yields in the 3rd channel
exo1 = 6.99
exo2 = 5.02
exo3 = 3.06

yield1 = 0.22
yield2 = 0.42
yield3 = 1. - yield1 - yield2

exo = yield1*exo1 + yield2*exo2 + yield3*exo3
ace_1d_exo_matrix, spindex.o2_p, spindex.e, exo, exomatrix

ace_1d_yieldmatrix, spindex.o2_p, spindex.e, spindex.o, 2.*yield1+yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.o2_p, spindex.e, spindex.o1d, yield2+2.*yield3, yieldmatrix

; O2+ + N(2D) -> NO+ + O / N+ + O2 ; Fox and Sung [2001]
exo = 0.
ace_1d_exo_matrix, spindex.o2_p, spindex.n2d, exo, exomatrix
yield1 = 1.8e-10/(1.8e-10 + 8.65e-11)
yield2 = 1. - yield1
ace_1d_yieldmatrix, spindex.o2_p, spindex.n2d, spindex.no_p, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.o2_p, spindex.n2d, spindex.o, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.o2_p, spindex.n2d, spindex.n_p, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.o2_p, spindex.n2d, spindex.o2, yield2, yieldmatrix

; O2+ + N(2P) -> N(4S) + O2+
; It is assumed here that the electronic energy goes into N(4S)
; Assuming value from Zipf 1980, which gives the transition wavelength
; of N(2P) -> N(4S) as 3466 Angstroms
exo = 3.57
ace_1d_exo_matrix, spindex.o2_p, spindex.n2p, exo, exomatrix

; O2+ + N(4S) -> NO+ + O
; Scott et al, JCP, 1998, Table II
exo  = 4.18
ace_1d_exo_matrix, spindex.o2_p, spindex.n4s, exo, exomatrix

; O2+ + NO -> NO+ + O2
; Midey et al JCP,(1999)
exo = 2.81
ace_1d_exo_matrix, spindex.o2_p, spindex.no, exo, exomatrix

; Neutral - Neutral Reactions
; ---------------------------
; 
; N(2P) + O2 -> NO + O
exo = 0.
ace_1d_exo_matrix, spindex.n2p, spindex.o2, exo, exomatrix

;N(2P) + O -> N(4S,2D) + O / NO+ + e
;See section 3.2.3 of Yonker [2013] Dissertation
; Value presently unknown
exo  = 0.000000000000000000000000000000000000000000
ace_1d_exo_matrix, spindex.n2p, spindex.o, exo, exomatrix
yield1 = 0.47 ; N(2D) channel
yield2 = 0.03 ; N(4S) channel
yield3 = 1. - yield1 - yield2
ace_1d_yieldmatrix, spindex.n2p, spindex.o, spindex.n2d, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n2p, spindex.o, spindex.n4s, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2p, spindex.o, spindex.no_p, yield3, yieldmatrix
ace_1d_yieldmatrix, spindex.n2p, spindex.o, spindex.e, yield3, yieldmatrix
ace_1d_yieldmatrix, spindex.n2p, spindex.o, spindex.o, yield1+yield2, yieldmatrix

; N(4S) + O2 -> NO + O
; Kennealy [1978], Duff[1994], Sharma [1998]
exo = 1.385
ace_1d_exo_matrix, spindex.n4s, spindex.o2, exo, exomatrix

; N(2D) + O2 -> NO + O(3P)
; Kennealy [1978], Rawlins [1989], Table 1. (1 kcal mol^-1 = 4.184 kJ mol^-1)
; This has to be treated correctly in NO chemiluminescence
exo = 3.76 ;
ace_1d_exo_matrix, spindex.n2d, spindex.o2, exo, exomatrix

; N(4S) + NO -> N2 + O 
; Sharma [1998] / Baulch [2005]
exo = 3.25
ace_1d_exo_matrix, spindex.n4s, spindex.no, exo, exomatrix

; N(2D) + NO -> N2 + O  
; Roble 1995
exo = 5.63
ace_1d_exo_matrix, spindex.n2d, spindex.no, exo, exomatrix

; N(2P) + O2 -> NO + O (3P,1D)
; Rawlins [1989], assuming only O(3P) channel. This has to be treated correctly in NO chemiluminescence
; Since vibrational yields of NO from this aren't known, we're good for now.
exo = 4.96
ace_1d_exo_matrix, spindex.n2p, spindex.o2, exo, exomatrix

; O + O + M -> O2 + M
; M is N2 or O2. Roble [1995]
exo = 5.12
ace_1d_exo_matrix, spindex.o, spindex.o, exo, exomatrix

; O + O2 + M -> O3 + M
; Herron, with t-dep
exo = 1.10
ace_1d_exo_matrix, spindex.o, spindex.o2, exo, exomatrix

; Quenching Reactions
; -------------------
; N(2D) + O -> N(4S) + O   
; Roble 1995
exo = 2.38
ace_1d_exo_matrix, spindex.n2d, spindex.o, exo, exomatrix

;N(2D) + N2 -> N(4S) + N2;
; Herron recommended
exo = 2.38
ace_1d_exo_matrix, spindex.n2d, spindex.n2, exo, exomatrix

; N(2D) + e -> N(4S) + e    
; Berrington and BUrke (1981)=Swaminathan, 1998 
; Richards [1986] showed that this energy goes into the electrons
exo = 2.38
ace_1d_exo_matrix, spindex.n2d, spindex.e, exo, exomatrix

; N(2P) + N2    will never be important 
exo = 3.57
ace_1d_exo_matrix, spindex.n2p, spindex.n2, exo, exomatrix

;N2P + NO
; Value presently unknown
exo = 0.0000000000000000000000000000000
ace_1d_exo_matrix, spindex.n2p, spindex.no, exo, exomatrix

; N2P + e -> N(4s,2d) + O (3P,1D)
; Herron, ignoring t-dep 
; Value presently unknown
exo = 0.0000000000000000000000000000000
ace_1d_exo_matrix, spindex.n2p, spindex.e, exo, exomatrix

; O(1D) + N2 -> O(3P) + N2
exo = 1.97
ace_1d_exo_matrix, spindex.o1d, spindex.n2, exo, exomatrix

; O(1D) + O2 -> O(3P) + O2
; This is assumed to be true, but a reference is needed for it
exo = 1.97
ace_1d_exo_matrix, spindex.o1d, spindex.o2, exo, exomatrix

; O(1D) + O -> O(3P) + O
; This is assumed to be true, but a reference is needed for it
exo = 1.97
ace_1d_exo_matrix, spindex.o1d, spindex.o, exo, exomatrix

; Auroral/Photoelectron absorption
;-----------------------------------
; N2 + e* -> 2N(4S)/2N(2D)/2N(2P) + e
yield1 = 0.5
yield2 = 0.276
yield3 = 1. - yield1 - yield2
ace_1d_yieldmatrix, spindex.n2, spindex.pe, spindex.n4s, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n2, spindex.pe, spindex.n2d, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2, spindex.pe, spindex.n2p, yield3, yieldmatrix
; The yields are assumed to be the same for auroral electrons
ace_1d_yieldmatrix, spindex.n2, spindex.ae, spindex.n4s, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n2, spindex.ae, spindex.n2d, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2, spindex.ae, spindex.n2p, yield3, yieldmatrix


; Photolysis of N2
yield1 = 0.5 ; N(4S)
yield2 = 0.275 ; N(2D)
yield3 = 1. - yield1 - yield2 ; N(2P)

ace_1d_yieldmatrix, spindex.n2, spindex.ph, spindex.n4s, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n2, spindex.ph, spindex.n2d, yield2, yieldmatrix
ace_1d_yieldmatrix, spindex.n2, spindex.ph, spindex.n2p, yield3, yieldmatrix
        
; Radiative Relaxation
;-----------------------------------
; The yields are written in terms of the relevant relaxation frequencies

; N(2D) -> N(4S) + hv   
exo = 2.38
ace_1d_exo_matrix, spindex.n2d, spindex.onebody, exo, exomatrix

; N(2P) -> N(2D)/N(4S) + hv
yield1 = 8.054e-2/(8.054e-2 + 5.308e-3)
yield2 = 5.308e-3/(8.054e-2 + 5.308e-3)
ace_1d_yieldmatrix, spindex.n2p, spindex.onebody, spindex.n2d, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.n2p, spindex.onebody, spindex.n4s, yield2, yieldmatrix

exo1 = 1.19
exo2 = 3.57
exo = yield1*exo1 + yield2*exo2

ace_1d_exo_matrix, spindex.n2p, spindex.onebody, exo, exomatrix

; O+(2P) -> O+(2D)/O+(4S) + hv
yield1 = 0.047/(0.171 + 0.047)
yield2 = 0.171/(0.171 + 0.047)
ace_1d_yieldmatrix, spindex.o2p_p, spindex.onebody, spindex.o_p, yield1, yieldmatrix
ace_1d_yieldmatrix, spindex.o2p_p, spindex.onebody, spindex.o2d_p, yield2, yieldmatrix

exo1 = 5.02*yield1
exo2 = 1.69*yield2
exo = yield1*exo1 + yield2*exo2
ace_1d_exo_matrix, spindex.o2p_p, spindex.onebody, exo, exomatrix

; O(1D) -> O(3P) + hv (630nm)
; Fennelly and Torr [1994]
exo = 1.97
ace_1d_exo_matrix, spindex.o1d, spindex.onebody, exo, exomatrix

