PRO ace_1d_chem_heating, heatmatrix, spindex, heatterms, model, i

    ev2ergs = 1.6022E-12
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ion-neutral reaction heating

q_in = fltarr(model.nlev)

; N+ + NO -> NO+ + N
    q_in = q_in + reform(heatmatrix[spindex.n_p, spindex.no,*]) ;*
    
; N+ + O -> O+ + N(4S)
    q_in = q_in + reform(heatmatrix[spindex.n_p, spindex.o,*]) ;*

; N+ + O2 -> NO+ + O                  
    q_in = q_in + reform(heatmatrix[spindex.n_p, spindex.o2,*]) ;*

; N2+ + N  -> N+ + N2       
    q_in = q_in + reform(heatmatrix[spindex.n2_p, spindex.n4s,*]) ;*

; N2+ + NO -> NO+ + N2
    q_in = q_in + reform(heatmatrix[spindex.n2_p,spindex.no,*]) ;*

; N2+ + O -> NO+ + N(2D)
    q_in = q_in + reform(heatmatrix[spindex.n2_p,spindex.o,*]) ;*

; N2+ + O2 -> O2+ + N2 ;
    q_in = q_in + reform(heatmatrix[spindex.n2_p,spindex.o2,*]) ;*

; O+ + N(2D) -> N+ + O 
    q_in = q_in + reform(heatmatrix[spindex.o_p, spindex.n2d,*]) ;*

; O+ + N2 -> NO+ + N(4S)                     
    q_in = q_in + reform(heatmatrix[spindex.o_p, spindex.n2,*]) ;*

; O+ + NO -> NO+ + O
    q_in = q_in + reform(heatmatrix[spindex.o_p, spindex.no,*]) ;*

; O+ + O2 -> O2+ + O
    q_in = q_in + reform(heatmatrix[spindex.o_p, spindex.o2,*]) ;*

; O+(2D)+ N2 -> N2+ + O
    q_in = q_in + reform(heatmatrix[spindex.o2d_p, spindex.n2,*])  ;*

; O+(2D) + O2 -> O2+ + O  
    q_in = q_in + reform(heatmatrix[spindex.o2d_p, spindex.o2,*]) ;*

; O+(2P) + N2 -> N2+ + O  
    q_in = q_in + reform(heatmatrix[spindex.o2p_p, spindex.n2,*]) ;*
 
; O2+ + N(4S) -> NO+ + O
    q_in = q_in + reform(heatmatrix[spindex.o2_p, spindex.n4s,*]) ;*

 ; O2+ + NO -> NO+ + O2
    q_in = q_in + reform(heatmatrix[spindex.o2_p, spindex.no,*]) ;*
  
    q_in = q_in * ev2ergs
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ion recombination heating
    q_ii = fltarr(model.nlev)

;N2+ + e -> N(2D) + N(4s)         
    q_ii = q_ii + reform(heatmatrix[spindex.n2_p, spindex.e,*]) ;*

; NO+ + e -> O + N(2D)/N(4S)
    q_ii = q_ii + reform(heatmatrix[spindex.no_p, spindex.e,*]) ;*

; O2+ + e -> O + O  
    q_ii = q_ii + reform(heatmatrix[spindex.o2_p, spindex.e,*]) ;*


    q_ii = q_ii * ev2ergs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Neutral-neutral reaction heating

    q_nn = fltarr(model.nlev)

; N(4S) + O2 -> NO + O
    q_nn = q_nn + reform(heatmatrix[spindex.n4s, spindex.o2,*]) ;*

; N(2D) + O2 -> NO + O     
    q_nn = q_nn + reform(heatmatrix[spindex.n2d, spindex.o2,*]) ;*

; N(4S) + NO -> N2 + O 
    q_nn = q_nn + reform(heatmatrix[spindex.n4s, spindex.no,*]) ;*

; N(2D) + NO -> N2 + O  
    q_nn = q_nn + reform(heatmatrix[spindex.n2d, spindex.no,*]) ;*

; N(2P) + O2 -> NO + O (3P,1D)
    q_nn = q_nn + reform(heatmatrix[spindex.n2p, spindex.o2,*]) ;*

; O + O + M -> O2 + M
    q_nn = q_nn + reform(heatmatrix[spindex.o, spindex.o, *]) ;*

; O + O2 + M -> O3 + M
    q_nn = q_nn + reform(heatmatrix[spindex.o, spindex.o2,*]) ;*
    
    q_nn = q_nn * ev2ergs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Quenching reaction heating

    q_quench = fltarr(model.nlev)

; O+(2D) + e -> O+(4S) + e + 3.31 ev ; Roble 95, table 1, reaction 26
    q_quench = q_quench + reform(heatmatrix[spindex.o2d_p, spindex.e, *]) ;*

; O+(2D) + O(3P) -> O+(4S) + O + 3.31 ev ; Roble 95, table 1, reaction 25
    q_quench = q_quench + reform(heatmatrix[spindex.o2d_p, spindex.o, *]) ;*
    
; O+(2P) + e -> O+(4S) / O+(2D)
    q_quench = q_quench + reform(heatmatrix[spindex.o2p_p, spindex.e, *]) ;*

; O+(2P) + O(3P) -> O+(4S) + O(3P) ;
    q_quench = q_quench + reform(heatmatrix[spindex.o2p_p, spindex.o, *]) ;*

; O2+ + N(2P) -> N(4S) + O2+
    q_quench = q_quench + reform(heatmatrix[spindex.o2_p, spindex.n2p, *]) ;*

; N(2D) + O -> N(4S) + O
    q_quench = q_quench + reform(heatmatrix[spindex.n2d, spindex.o, *]) ;*
    
;N(2D) + N2 -> N(4S) + N2;
    q_quench = q_quench + reform(heatmatrix[spindex.n2d, spindex.n2, *]) ;*
    
; N(2D) + e -> N(4S) + e  ; Richards [1986] showed that this energy goes into the electrons
; This has been accounted for in ace_1d_eionheating
    q_quench = q_quench + 0.*reform(heatmatrix[spindex.n2d, spindex.e, *]) ;*
    
; O(1D) + N2 -> O + N2
    q_quench = q_quench + reform(heatmatrix[spindex.o1d, spindex.n2,*]) ;*
    
; O(1D) + O2 -> O + O2
    q_quench = q_quench + reform(heatmatrix[spindex.o1d, spindex.o2,*]) ;*
    
; O(1D) + O -> O + O
    q_quench = q_quench + reform(heatmatrix[spindex.o1d, spindex.o,*])
    
    q_quench = q_quench * ev2ergs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Airglow processes that remove energy

	q_airglow = fltarr(model.nlev)

; N(2D) -> N(4S) + hv 
	q_airglow = q_airglow + reform(heatmatrix[spindex.n2d, spindex.onebody, *]) ;*
	
; N(2P) -> N(2D)/N(4S) + hv
	q_airglow = q_airglow + reform(heatmatrix[spindex.n2p, spindex.onebody, *]) ;*
	
; O+(2P) -> O+(4S)/O+(2D) + hv
	q_airglow = q_airglow + reform(heatmatrix[spindex.o2p_p, spindex.onebody, *]) ;*
	
; O(1D) -> O(3P) + hv (630nm)
	q_airglow = q_airglow + reform(heatmatrix[spindex.o1d, spindex.onebody, *]) ;*
	
	q_airglow = q_airglow * ev2ergs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    heatterms.q_ionneut  = q_in
    heatterms.q_ionrec   = q_ii
    heatterms.q_neutneut = q_nn
    heatterms.q_quench   = q_quench
    heatterms.q_airglow = q_airglow
    heatterms.q_chem = q_ii + q_in + q_nn + q_quench - q_airglow
    heatterms.q_total = heatterms.q_chem ; q_total is re-initalized here, and q_chem is the first term in the expression
    
END
