; Parameters describing the model run are defined here

; Model constants
model  = {  del_time : inputs.timestep      ,$ ; Time step, in seconds
            nlev     : 57.                  ,$ ; 57 vertical grid points
            ztop     : 7.                   ,$ ; +7 is the top-most vertical grid point (pressure level)
            zbot     :-7.                   ,$ ; -7 is the bottom-most vertical grid point
            p0       : 5.e-4                ,$ ; Model reference pressure in dyne/cm^2, corresponding to the pressure at the grid point z = 0          
            zp       : findgen(57)/4 - 7.   ,$ ; The vertical grid has 4 grid points per scale height. 
            zbound   : 97e5                 ,$ ; The lower pressure boundary of the model is assumed to be at 97 km
            dz       : 0.25                 ,$ ; dz = 0.25

	         nday_hours      : 24, $
	         nhour_minutes   : 60, $
	         nminute_seconds : 60  $
         }            

; Physical constants and conversion factors
pconst = {  boltz : 1.38066e-16,$       ; boltzmann constant  in ergs K^-1        
            h     : 6.626e-27  ,$       ; in ergs s
            c     : 2.9979e10  ,$       ; in cm s^-1
            avo   : 6.023e23   ,$       ; avogadro's number           
            gask  : 8.314e7    ,$       ; gas constant                
            grav  : fltarr(model.nlev),$  ; acceleration due to gravity varies as a function of height, in cm s^-2
            gravref: 980.      ,$       ; Reference value of gravity
            re    : 6371.e5    ,$       ; Earth radius in cm
            cgm   : 0.5        ,$       
            eff   : 0.33       ,$       ; Solar EUV heating efficiency, from TIE-GCM
            ev2erg: 1.602E-12  ,$		; eV to ergs
            amu   : 1.6605e-24 ,$
            q_e   : 1.602e-19  ,$       ; Electron charge
            m_e   : 9.109e-28   $		; Electron mass
         }


; Atomic/Molecular molar masses   
mass =  {   o2    : 32. , $
            o     : 16. , $
            n2    : 28. , $
            ar    : 40. , $ 
            na    : 23. , $
            he    : 4.  , $
            no    : 30. , $ 
            n4s   : 14. , $
            n2d   : 14. , $
            ch4   : 16. , $
            h2    : 2.  , $
            co    : 28. , $
            co2   : 44. , $
            h2o   : 18. , $ 
            hox   : 1.  , $
            h     : 1.    $
        }
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Chemical species present in the model;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;         

;; This value *must* be the same as in ace_1d_update_crh
nsp=23 ; Number of species defined in the model -- though not all of them are presently used

spindex={$
  nsp:  nsp, $
  O:      0, $
  O2:     1, $
  N2:     2, $
  NO:     3, $
  N4S:    4, $
  N2D:    5, $
  N2P:    6, $
  N2A:    7, $
  O1D:    8, $
  O3P:    9, $
  O_P:   10, $
  O2_P:  11, $
  N2_P:  12, $
  NO_P:  13, $
  N_P :  14, $
  O2D_P: 15, $
  e:     16, $
  O3:    17, $
  O2P_P:    18, $    ; Addition by KV 6/15/16
  onebody: nsp-4, $  ; For reactions involving radiative relaxation
  ae:      nsp-3, $
  pe:      nsp-2, $
  ph:      nsp-1}

spname = { $
  O:     'O'        	 , $
  O2:    'O!d2!n'   	 , $
  N2:    'N!d2!n'   	 , $
  NO:    'NO'       	 , $
  N4S:   'N(!u4!nS)'	 , $
  N2D:   'N(!u2!nD)'	 , $
  N2P:   'N(!u2!nP)'	 , $
  N2A:   'N!d2!nA)' 	 , $
  O1D:   'O!u1!nD'  	 , $
  O3P:   'O!u3!nP'  	 , $
  O_P:   'O!u+!n'   	 , $
  O2_P:  'O!d2!u+!n'	 , $
  N2_P:  'N!d2!u+!n'	 , $
  NO_P:  'NO!u+!n'  	 , $
  N_P:   'N!u+!n'   	 , $
  O2D_P: 'O(!u2!nD!u+!n)', $
  e:     'e!u-!n'   	 , $
  O3:     'O!d3!n'       , $
  O2P_P: 'O(!u2!nP!u+!n)', $
  onebody:'1Body'    	 , $  
  ae:     'Ae'       	 , $
  pe:     'Pe'       	 , $
  ph:     'Ph'         }
