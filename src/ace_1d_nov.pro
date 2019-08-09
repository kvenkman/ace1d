FUNCTION ace_1d_nov, zmaj, zminor, zion, termsno, model, pconst

av1 = fltarr(11, 11)
av2 = fltarr(11, 11)
knet  = fltarr(11, 11)
knetp = fltarr(11, 11)
lossmatrix = fltarr(11, 11)
nov = fltarr(model.nlev)

n4s_prod = termsno.no_n4s_prod; NO production rate from N(4S)
n2d_prod = termsno.no_n2d_prod; NO production rate from N(2D)
no_loss = termsno.no_loss

; Defining the array av1:
      av1[*, 0]  = [0.,12.396,0.,0.,0.,0.,0.,0.,0.,0.,0. ]
      av1[*, 1]  = [0.,0.,23.405,0.,0.,0.,0.,0.,0.,0.,0. ]
      av1[*, 2]  = [0.,0.,0.,33.122,0.,0.,0.,0.,0.,0.,0. ]
      av1[*, 3]  = [0.,0.,0.,0.,41.625,0.,0.,0.,0.,0.,0. ] 
      av1[*, 4]  = [0.,0.,0.,0.,0.,49.006,0.,0.,0.,0.,0. ]
      av1[*, 5]  = [0.,0.,0.,0.,0.,0.,55.344,0.,0.,0.,0. ]
      av1[*, 6]  = [0.,0.,0.,0.,0.,0.,0.,60.723,0.,0.,0. ]
      av1[*, 7]  = [0.,0.,0.,0.,0.,0.,0.,0.,65.212,0.,0. ]
      av1[*, 8]  = [0.,0.,0.,0.,0.,0.,0.,0.,0.,68.876,0. ]
      av1[*, 9]  = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,71.789 ]
      av1[*, 10] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.     ]

; Defining the array av2:
      av2[*, 0]  = [0.,0.,0.786,0.,0.,0.,0.,0.,0.,0. ,0. ]
      av2[*, 1]  = [0.,0.,0.,2.281,0.,0.,0.,0.,0.,0. ,0. ]
      av2[*, 2]  = [0.,0.,0.,0.,4.400,0.,0.,0.,0.,0. ,0. ] 
      av2[*, 3]  = [0.,0.,0.,0.,0.,7.052,0.,0.,0.,0. ,0. ] 
      av2[*, 4]  = [0.,0.,0.,0.,0.,0.,10.141,0.,0.,0.,0. ] 
      av2[*, 5]  = [0.,0.,0.,0.,0.,0.,0.,13.574,0.,0.,0. ]  
      av2[*, 6]  = [0.,0.,0.,0.,0.,0.,0.,0.,17.257,0.,0. ] 
      av2[*, 7]  = [0.,0.,0.,0.,0.,0.,0.,0.,0.,21.104,0. ] 
      av2[*, 8]  = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,25.036 ] 
      av2[*, 9]  = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.     ]
      av2[*, 10] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.     ]

; Defining the array knet:
;	  Collisional relaxation rates obtained from caridade 2008
;     rates for relaxation from v=5,7,9 not available,interpolated from 
;     k4,k6,k8,k10

;@298 K 
; For v = 1, using the Hwang rate at 298 K, and the Caridade rate at 1500K
; Changed the v=1 rate to the Caridade rate (KV 11/01/17)
      knet[*, 0] =[0.,21.24,15.58,11.25,8.33,0.,6.08,0.,5.76,0.,4.66]  
      knet[*, 1] =[0.,0.   ,9.20 ,8.76 ,7.95,0.,6.90,0.,4.84,0.,4.74]
      knet[*, 2] =[0.,0.   ,0.   ,6.26 ,6.79,0.,5.41,0.,4.56,0.,4.57]
      knet[*, 3] =[0.,0.   ,0.   ,0.   ,4.68,0.,5.01,0.,4.07,0.,4.21]
      knet[*, 4] =[0.,0.   ,0.   ,0.   ,0.  ,0.,4.16,0.,4.01,0.,3.83]
      knet[*, 5] =[0.,0.   ,0.   ,0.   ,0.  ,0.,3.05,0.,3.20,0.,2.52]
      knet[*, 6] =[0.,0.   ,0.   ,0.   ,0.  ,0.,0.  ,0.,3.18,0.,2.77]
      knet[*, 7] =[0.,0.   ,0.   ,0.   ,0.  ,0.,0.  ,0.,2.47,0.,2.58]
      knet[*, 8] =[0.,0.   ,0.   ,0.   ,0.  ,0.,0.  ,0.,0.  ,0.,2.42]
      knet[*, 9] =[0.,0.   ,0.   ,0.   ,0.  ,0.,0.  ,0.,0.  ,0.,2.19]
      knet[*, 10]=[0.,0.   ,0.   ,0.   ,0.  ,0.,0.  ,0.,0.  ,0.,0.  ]
;@1500 K       
      knetp[*, 0] =[0.,20.71,15.14,0. ,8.98 ,0.,6.98,0.,6.30,0.,5.48]  
      knetp[*, 1] =[0.,0.   ,11.91,0. ,8.95 ,0.,7.31,0.,5.71,0.,5.25]
      knetp[*, 2] =[0.,0.   ,0.   ,0. ,7.53, 0.,6.16,0.,5.84,0.,4.94]
      knetp[*, 3] =[0.,0.   ,0.   ,0.  ,6.96,0.,6.13,0.,5.04,0.,4.48]
      knetp[*, 4] =[0.,0.   ,0.   ,0.  ,0.  ,0.,4.97,0.,4.75,0.,4.27]
      knetp[*, 5] =[0.,0.   ,0.   ,0.  ,0.  ,0.,5.14,0.,4.59,0.,3.76]
      knetp[*, 6] =[0.,0.   ,0.   ,0.  ,0.  ,0.,0.  ,0.,3.89,0.,3.73]
      knetp[*, 7] =[0.,0.   ,0.   ,0.  ,0.  ,0.,0.  ,0.,3.62,0.,3.62]
      knetp[*, 8] =[0.,0.   ,0.   ,0.  ,0.  ,0.,0.  ,0.,0.  ,0.,3.06]
      knetp[*, 9] =[0.,0.   ,0.   ,0.  ,0.  ,0.,0.  ,0.,0.  ,0.,2.82]
      knetp[*, 10]=[0.,0.   ,0.   ,0.  ,0.  ,0.,0.  ,0.,0.  ,0.,0.  ]
  
	  knet[0,0]  = 1.e-10
 	  knetp[0,0] = 1.e-10
; Filling in for v=5,7,9 by interpolation
      knet[5,*]  = (knet[4,*] + knet[6, *])/2.
      knet[7,*]  = (knet[6,*] + knet[8, *])/2.
      knet[9,*]  = (knet[8,*] + knet[10, *])/2.

; Filling in for v=3,5,7,9 by interpolation
      knetp[3, *]  = (knetp[4, *] + knetp[2, *])/2.
      knetp[5, *]  = (knetp[4, *] + knetp[6, *])/2.
      knetp[7, *]  = (knetp[6, *] + knetp[8, *])/2.
      knetp[9, *]  = (knetp[8, *] + knetp[10, *])/2.

; Collisional relaxation from & into the same level should be zero
       knetp[3,3] = 0.
       knetp[5,5] = 0.
       knetp[7,7] = 0.
       knetp[9,9] = 0.
       
       knet[5,5] = 0.
       knet[7,7] = 0.
       knet[9,9] = 0.
       
; For interpolating knet to neutral temperature used at current step       
; Assuming the form y=mx+c, interp gives the term 'c'
      interp  = knetp - 1500.*((knetp - knet)/1202.)

; N(4S) yields estimated from Sultanov & Balakrishnan (2006), Fig. 5
      g_4s = [0.106307,0.116395,0.117199,0.147256,0.168128, $
       0.127784,0.101884,0.0682910,0.0292211,0.0122811,0.00525504]
; N(2D) yields estimated from Miquel (2003)
       g_2d = [0.0552162,0.0520401,0.0508185,0.0408014,0.0522844, $
       0.0752504,0.111654,0.135353,0.144149,0.121915,0.160518]
    
; Loop over altitude
    
    FOR i = 0, model.nlev - 1 DO BEGIN     
      knetf = interp + zmaj.tn[i]*((knetp - knet)/1202.)
      knetf = knetf*1e-12

      prodmatrix = g_4s*n4s_prod[i] + g_2d*n2d_prod[i]
      lossmatrix=-1.*(av1+av2+(knetf*zmaj.oden[i])) ; "Negative loss"
      
      FOR j = 0, 10 DO BEGIN
        lossmatrix[j, j] = lossmatrix[j, j] + total(av1[j, *]) + $
                           total(av2[j, *]) + total(knetf[j, *])*zmaj.oden[i] + $
                           no_loss[i]
      ENDFOR
             
      nov_levels = la_linear_equation(lossmatrix, prodmatrix, status=status)
      
      nov[i] = (pconst.h*pconst.c/5.3e-4)*( $
                nov_levels[1]*av1[1, 0] + nov_levels[2]*av1[2, 1] + nov_levels[3]*av1[3, 2] + $
                nov_levels[4]*av1[4, 3] + nov_levels[5]*av1[5, 4] + nov_levels[6]*av1[6, 5] + $
                nov_levels[7]*av1[7, 6] + nov_levels[8]*av1[8, 7] + nov_levels[9]*av1[9, 8] + $
                nov_levels[10]*av1[10, 9])+ $
               (pconst.h*pconst.c/2.7e-4)*( $
                nov_levels[2]*av2[2, 0] + nov_levels[3]*av2[3, 1] + nov_levels[4]*av2[4, 2] + $
                nov_levels[5]*av2[5, 3] + nov_levels[6]*av2[6, 4] + nov_levels[7]*av2[7, 5] + $
                nov_levels[8]*av2[8, 6] + nov_levels[9]*av2[9, 7] + nov_levels[10]*av2[10, 8])
      
      IF(status > 0) then print, "Uh oh Uh oh Uh oh Uh oh Uh oh Uh oh Uh oh. Something broke in novcalc"        
    ENDFOR
      
    ;STOP 
    return, nov; in ergs cm^-3
END
