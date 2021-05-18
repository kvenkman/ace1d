; Setting up the variables

             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;            
;;;;    General Purpose structures    ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zmaj = {    o      : fltarr(model.nlev)   ,$
            o2     : fltarr(model.nlev)   ,$
            n2     : fltarr(model.nlev)   ,$
            oden   : fltarr(model.nlev)   ,$
            o2den  : fltarr(model.nlev)   ,$
            n2den  : fltarr(model.nlev)   ,$
            tn     : fltarr(model.nlev)   ,$
            z      : fltarr(model.nlev)   ,$ ; Altitude of pressure levels
            zz     : fltarr(model.nlev)   ,$ ; Altitude in km
            barm   : fltarr(model.nlev)   ,$ ; mean mass
            cp     : fltarr(model.nlev)   ,$
            kt     : fltarr(model.nlev)    $
       }

; minor constituents in number density units           
zminor = {  no      : fltarr(model.nlev)   ,$  ; NO with diffusion
            no_mmr  : fltarr(model.nlev)   ,$
            n4s_mmr : fltarr(model.nlev)   ,$
            no_pce  : fltarr(model.nlev)   ,$  ; NO in PCE
            co2     : fltarr(model.nlev)   ,$
            n4s     : fltarr(model.nlev)   ,$
            n4s_pce : fltarr(model.nlev)   ,$
            n2d     : fltarr(model.nlev)   ,$ ;N(2D)
            n2p     : fltarr(model.nlev)   ,$ ;N(2P)
            n2a     : fltarr(8, model.nlev),$ ;N2(A) for v = 0-7
            o1d     : fltarr(model.nlev)   ,$
            ; These species aren't currently used in the model
            oh      : fltarr(model.nlev)   ,$
            o3      : fltarr(model.nlev)   ,$
            no2     : fltarr(model.nlev)   ,$
            h2o     : fltarr(model.nlev)   ,$
            h       : fltarr(model.nlev)   ,$
            ho2     : fltarr(model.nlev)   ,$
            ox      : fltarr(model.nlev)   ,$
            nox     : fltarr(model.nlev)   ,$
            hox     : fltarr(model.nlev)   ,$
            noxno   : fltarr(model.nlev)   ,$
            no3     : fltarr(model.nlev)   } ; array of minor species

zion = {    o_p    : fltarr(model.nlev)   ,$
            n_p    : fltarr(model.nlev)   ,$
            no_p    : fltarr(model.nlev)  ,$
            o_p_pce : fltarr(model.nlev)  ,$
            n_p_pce : fltarr(model.nlev)  ,$
            no_p_pce : fltarr(model.nlev) ,$
            o2d_p   : fltarr(model.nlev)  ,$            
            o2p_p   : fltarr(model.nlev)  ,$  ; Addition by KV 6/15/16
            n2_p    : fltarr(model.nlev)  ,$  ; N2+
            o2_p    : fltarr(model.nlev)  ,$            
            e       : fltarr(model.nlev)  ,$
            te      : fltarr(model.nlev)  ,$
            ti      : fltarr(model.nlev)   } ; array of ions
            
            
model_time = { day   : fltarr(1) ,$
               time  : fltarr(1) ,$
               year  : fltarr(1)  $
              }
              
model_sun = { f107d : 0.0 ,$ ; Daily F10.7cm flux
              f107a : 0.0 ,$ ; 81-day average F10.7cm flux
              kp    : fltarr(8) ,$ 
              ap    : 0.0,$ 
              sza   : 60.0,$ ; Solar Zenith angle
              ecc   : 1.0 $ ; eccentricity factor due to orbit
            }

; Photoelectron processes
; PE = Photoelectron
; E.g. ped_o2 = photoelectron dissociation of O2
; aei_n2 = ionization of N2 by auroral electrons
; NOTE: Auroral processes defined but not used presently 04/10/2019

zei     =   {$
            pei_n2:fltarr(model.nlev) 		,$
            pei_o2:fltarr(model.nlev) 		,$
            pei_o:fltarr(model.nlev)  		,$
            pei_o_op2p : fltarr(model.nlev) ,$
            pei_o_op2d : fltarr(model.nlev) ,$
            pei_o_op4s : fltarr(model.nlev) ,$  
            di_o2 : fltarr(model.nlev)	    ,$
            di_n2 : fltarr(model.nlev)		,$
            ped_o2:fltarr(model.nlev) 		,$                     
            ped_n2:fltarr(model.nlev) 		,$
            aei_n2:fltarr(model.nlev) 		,$
            aei_o2:fltarr(model.nlev) 		,$
            aei_o:fltarr(model.nlev)  		,$
            aed_n2:fltarr(model.nlev)  		 $
            }

; Solar EUV driven ionization processes            
zpid    =    {$
             i_o:fltarr(model.nlev)  ,$
             i_o_op2p : fltarr(model.nlev), $
             i_o_op2d : fltarr(model.nlev), $
             i_o_op4s : fltarr(model.nlev), $            
             i_n2:fltarr(model.nlev) ,$
             i_o2:fltarr(model.nlev) ,$
             i_no:fltarr(model.nlev) ,$
             i_n4s:fltarr(model.nlev)  ,$
             j_n2:fltarr(model.nlev) ,$
             j_o2:fltarr(model.nlev) ,$
             j_o2_o1d:fltarr(model.nlev) ,$
             j_o2wvln:fltarr(model.nlev,solspec.n_wave), $            
             j_no:fltarr(model.nlev) ,$
             di_n2:fltarr(model.nlev),$
             di_o2:fltarr(model.nlev) $
             }
             
zcol    =   {$
            o   :fltarr(model.nlev) ,$ ; vertical column density
            o2  :fltarr(model.nlev) ,$ ; vertical column density
            n2  :fltarr(model.nlev) ,$ ; vertical column density
            so   :fltarr(model.nlev),$ ; slant integral for current zenith angle
            so2  :fltarr(model.nlev),$ ; slant integral for current zenith angle
            sn2  :fltarr(model.nlev),$ ; slant integral for current zenith angle
            w    :fltarr(model.nlev) $ ; altitude difference between pressure levels
            }           


; Structure containing energy sources
edep = {$
            q          : fltarr(model.nlev), $
            n2_abs     : fltarr(model.nlev)  $
       }

; Structure containing energy loss sources
eloss = {$
            explicit          : fltarr(model.nlev), $
            implicit          : fltarr(model.nlev)  $
        }
        
; Molecular diffusion coefficients
k_m       = {$

            no  : fltarr(model.nlev), $
            n4s : fltarr(model.nlev), $
            o   : fltarr(model.nlev), $
            o2  : fltarr(model.nlev), $
            n2   : fltarr(model.nlev) $
            
            }
            
flux = { $
            no  : fltarr(model.nlev), $
            n4s : fltarr(model.nlev), $
            o   : fltarr(model.nlev), $
            o2  : fltarr(model.nlev), $
            n2  : fltarr(model.nlev)  $
}

k_e = { k : 5e-6*exp(-7-model.zp)}

coolterms = { $
             noc:fltarr(model.nlev), $
             nov:fltarr(model.nlev), $
             co2_cool:fltarr(model.nlev), $
             o3p_cool:fltarr(model.nlev), $
             net_cool:fltarr(model.nlev) $
            }
            
heatterms = { $
            q_neutneut : fltarr(model.nlev), $
            q_ionrec   : fltarr(model.nlev), $
            q_ionneut  : fltarr(model.nlev), $
            q_quench   : fltarr(model.nlev), $
            q_airglow  : fltarr(model.nlev), $
            q_chem     : fltarr(model.nlev), $ ; this is the sum of the above four terms
            q_euv      : fltarr(model.nlev), $
            q_thermale : fltarr(model.nlev), $
            q_thermali : fltarr(model.nlev), $
            q_srb      : fltarr(model.nlev), $
            q_src      : fltarr(model.nlev), $
            q_joule    : fltarr(model.nlev), $
            q_aurora   : fltarr(model.nlev), $
            q_total    : fltarr(model.nlev), $
            s_ped      : fltarr(model.nlev)  $
            }
cooltermse = { $
            n2rot  : fltarr(model.nlev), $
            o2rot  : fltarr(model.nlev), $
            co2rot : fltarr(model.nlev), $
            h2orot : fltarr(model.nlev), $
            n2vib  : fltarr(model.nlev), $
            o2vib  : fltarr(model.nlev), $
            ofine  : fltarr(model.nlev), $
            o1d    : fltarr(model.nlev), $
            elasticn2 : fltarr(model.nlev), $
            elastico2 : fltarr(model.nlev), $
            elastico  : fltarr(model.nlev), $
            ei        : fltarr(model.nlev), $
            explicit :  fltarr(model.nlev), $
            neutrals_implicit :  fltarr(model.nlev), $
            ions_implicit :  fltarr(model.nlev), $
            total  : fltarr(model.nlev)  $
            }
            
heattermse = { $
            q_pe0_55   : fltarr(model.nlev), $
            q_pe55_105 : fltarr(model.nlev), $
            q_quench : fltarr(model.nlev), $
            q_total       : fltarr(model.nlev) $
            }

cooltermsi = { $
             neut:fltarr(model.nlev) $
            }
            
heattermsi = { $
            ei   : fltarr(model.nlev), $
            q_pe : fltarr(model.nlev), $
            q_ee   : fltarr(model.nlev) $
            }

o_o2_terms = {$
            s11 : fltarr(model.nlev), $
            s12 : fltarr(model.nlev), $
            s21 : fltarr(model.nlev), $
            s22 : fltarr(model.nlev), $
            s10 : fltarr(model.nlev), $
            s20 : fltarr(model.nlev)  $
            }
            
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Chemistry related structures;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   yieldmatrix  = fltarr(spindex.nsp,spindex.nsp,spindex.nsp)

   exomatrix    = fltarr(spindex.nsp,spindex.nsp)

   coeffmatrix = fltarr(spindex.nsp,spindex.nsp,model.nlev)
   ratematrix  = fltarr(spindex.nsp,spindex.nsp,model.nlev)
   heatmatrix  = fltarr(spindex.nsp,spindex.nsp,model.nlev)
   
   p_ave = fltarr(4, 8)
