PRO ace_1d_ratecoeff, zmaj, zion, zminor, model, chemk, chema,  $
                      chemg, chemq, chemr, chemf, ratek, ratea, $
                      rateg, rateq, rater, p_ave
;REACTION COEFFICIENTS procedure, taken from 1D NO model as of January, 2016
; Scott and Karthik

; K                   : ion-neutral     
; G (gamma)           : neutral-neutral      
; A (alpha)           : ion recombination    
; Q                   : quenching      
; R                   : radiative relaxation   
; F                   : branching ratio     
;** indices in each variabile name (e. g. Q50) refer to:
;***  O................................ 0
;***  O2................................1
;***  N2................................2
;***  NO................................3
;***  N(4S).............................4
;***  N(2D).............................5
;***  e.................................6
;***  O+(2D)............................7
;***  N2(A(.............................8
;***  N(2P).............................9
;O(1D).................................10
;O(3P) (treated as minor species)......11
;O3....................................12
; 
;
;  Ion Neutral Reactions
;  ---------------------
;  
;  O+(2D) + O(3P) - O+(4S) + O(3P)
chemk.k_00 = 5.0e-11
ratek.k_00 = chemk.k_00 * zion.o2d_p * zmaj.oden
; Stephan et al 2003; see Justin's thesis

;O+ + O2 -> O2+ + O
k_01_label = 'O!u+!n + O!d2!n -> O!d2!u+!n + O'
chemk.K_01 = 1.7e-11*(300./zmaj.tn)^(.77) +8.54e-11*exp(-28.8/4.184/23.06/1.38e-23/zmaj.tn/(1.0/1.6e-19));
;from Hierl et al, JCP, 1997.  The constants--unelegantly-- convert from kj/mol (as listed by Hierl) to kcal/mol to  ev to joules = kT)
ratek.k_01= chemk.k_01 * zion.o_p * zmaj.o2den

; O+ + N2 -> NO+ + N(4S)                 
k_02_label = 'O!u+!n + N!d2!n -> NO+ + N(!u4!nS)'    
chemk.K_02 = 1.20E-12*(300./zmaj.tn)^.45
;from Hierl et al (1997), courtesy of Fox and Sung, JGR, 2001
ratek.k_02 = chemk.k_02 * zion.o_p * zmaj.n2den

; O+ + NO -> NO+ + O
chemk.k_03 =  5.01e-13*(300./zmaj.tn)^(1.68) + 4.02e-12*exp(-7.5/4.184/23.06/1.38e-23/zmaj.tn/(1/1.6e-19)) 
; new reaction added 9/2008 by jy.  reference:  Dotan and Viggiano, JCP, 1999
ratek.k_03 = chemk.k_03 * zion.o_p * zminor.no
           
; O2+ + NO -> NO+ + O 
chemk.k_13 = 4.5e-10
; Midey etal JCP,(1999)--independent of temperature
ratek.k_13 = chemk.k_13 * zion.o2_p * zminor.no
       
; O2+ + N(4S) -> NO+ + O
chemk.K_14 = 1.33e-10
;Average of measurements cited by  Scott et al, JCP, 1998
ratek.k_14 = chemk.k_14 * zion.o2_p * zminor.n_4s       
;STOP
 ; N2+ + O -> NO+ + N(2D)       
 chemk.K_20 = 1.4E-10*(zmaj.tn/300.)^(-0.44)
 ;         -> O+ + N2 same as FOx and Sung (2001) when BR(N2d) = .95
ratek.k_20 = chemk.k_20 * zion.n2_p * zmaj.oden

; N2+ + O2 -> O2+ + N2 (orig)    
chemk.K_21 = 5.1E-11*(300./zmaj.tn)^(1.16)
; see Fox and Sung (2001)
ratek.k_21 = chemk.k_21 * zion.n2_p * zmaj.o2den
 
;  N2+ + NO -> NO+ + N2       
chemk.k_23 = 7.5e-9*(zmaj.tn)^(-.52)
; new reaction added 9/2008 by jy.  reference Midey, Miller, Viggiano, JCP, 2004.
ratek.k_23 = chemk.k_23 * zion.n2_p * zminor.no

; O+(2D) + O2 -> O2+ + O  
chemk.k_71 = 7.e-10
;from Johnsen and Biondi, JCP, (1980)
ratek.k_71 = chemk.k_71 * zion.o2d_p * zmaj.o2den
       
; O+(2D)+N2 -> N2+ + O
chemk.K_72 = 5.7e-10*exp(-400./zmaj.tn)
;see Fox and Sung (2001)
ratek.k_72 = chemk.k_72 * zion.o2d_p * zmaj.n2den

; N+ + O -> O+ + N(4S)
chemk.K_40 = 4.5e-12
;from Anicich, JPL, 2003        
ratek.k_40 = chemk.k_40 * zion.n_p * zmaj.oden
     
; N+ + O2 -> NO+ + O                  
chemk.K_41 = 5.5e-10 
;Midey et al, JCP 2006. ; See branching ratios below under F_5 and F_6.
ratek.k_41 = chemk.k_41 * zion.n_p * zmaj.o2den
           
; N + N2+ -> N+ + N2       
chemk.k_42 = 1.e-11
; Fox and SUng (2001)
ratek.k_42 = chemk.k_42 * zminor.n_4s * zion.n2_p
                    
; N+ + NO -> NO+ + N       
chemk.k_43 = 6.5e-9*(zmaj.tn)^(-.44)
;new reaction added 9/2008 by jy.  reference Midey, Miller, Viggiano, JCP, 2004
ratek.k_43 = chemk.k_43 * zion.n_p * zminor.no

;
; Ion Recombination
; ---------------------------
; 
; where is O+ recombination A0?

; O2+ + e -> O + O  
chema.a_1 = 1.95E-7*(zmaj.tn/300.)^(-.7)
;Sheehan and St.Maurice, jgr (2004)
ratea.a_1 = chema.a_1 * zion.o2_p * zion.elec       
  
;N2+ + e -> N(2D) + N(4s)         
chema.a_2 = 2.2E-7*(zmaj.tn/300.)^(-.39)
;Sheehan and St.maurice, jgr,(2004)
ratea.a_2 = chema.a_2 * zion.n2_p * zion.elec       

;(NO+ + e -> O + N(4S) 
; NO+ + e -> O + N(2D)      
chema.a_3 = 3.5E-7*(zmaj.tn/300.)^(-.69)
;Sheehan and ST.Maurice(2004)
ratea.a_3 = chema.a_3 * zion.no_p * zion.elec

; N2+ + e -> N2(A), 8 vibrational levels?
chema.a_8 = [0.422106,0.377770,0.347252,0.323350,0.311022,0.305038,0.299031,0.422961]    


;
; Neutral - Neutral Reactions
; ---------------------------
; 
; N(4S) + O2 -> NO + O
chemg.g_14 = 1.5E-11*exp(-3600./zmaj.tn)
; JPL? Need to check Justin's thesis to get the orgigin of this rate
rateg.g_14 = chemg.g_14 * zminor.n_4s * zmaj.o2den

; N(2D) + O2 -> NO + O     
chemg.g_15 = 6.2e-12*(zmaj.tn/300.)
;Duff (2003) (t3)
rateg.g_15 = chemg.g_15 * zminor.n_2d * zmaj.o2den
       
; N(4S) + NO -> N2 + O 
chemg.g_34 = 2.1e-11*exp(100./(zmaj.tn))   
;this is the jpl 2006 recommendation
rateg.g_34 = chemg.g_34 * zminor.n_4s * zminor.no

; N(2D) + NO -> N2 + O  
chemg.g_35 = 6.E-11;  (t3)          
;Herron (1999)
rateg.g_35 = chemg.g_35 * zminor.n_2d * zminor.no

;N2P + O -> N(4S,2D) + O (3P,1D);  BR in F_16
chemg.g_90=2.7e-11
;???Heron???
rateg.g_90 = chemg.g_90 * zminor.n_2p * zmaj.oden

;N2P + O2 -> NO + O (3P,1D)
chemg.g_91 = 3.1e-12*exp(-60./zmaj.tn)
;Herron, with t-dep
rateg.g_91 = chemg.g_91 * zminor.n_2p * zmaj.o2den

;
; Quenching Reactions
; -------------------

; N(2D) + O -> N(4S) + O   
chemq.q_50=1.65E-12*exp(-260./zmaj.tn)
; Fell with Davenport TD
;JY 7/7--Herron(1999) somehow omits self-described 'definitive' lab work on this reaction, 
;i.e. Fell and Steinfeld, JCP, 1990.  Herron fit the temperature dependence found by Davenport et al, JGR,81,1976 
;to an average of the reported room temperature rate coefficients.  The same is done here:  the temperature dependence of Davenport is fit so that at TEMP=298 K, Q_50=6.9e-13, 
; which is the Fell and Steinfeld result. 
rateq.q_50 = chemq.q_50 * zminor.n_2d * zmaj.oden

;N(2D) + N2 -> N(4S) + N2;
chemq.Q_52=1.74e-14; Herron recommended
rateq.q_52 = chemq.q_52 * zminor.n_2d * zmaj.n2den

;N(2D) + e -> N(4S) + e    
chemq.q_56=3.86E-10*(zmaj.tn/300.)^.81 ;Berrington and BUrke (1981)=Swaminathan, 1998 
rateq.q_56 = chemq.q_56 * zminor.n_2d * zion.elec

;Total deactivation rate of n2a(v=3-7) by o.  Branching ratio of
;.57 is the amount into the NO + N(4S,2D)channel (Thomas and Kaufman, JPC, 1997).  Temp
;dependence from Hill JGR,1999.
chemq.q_80=[2.80000e-11,3.40000e-11,3.34000e-11,3.48000e-11, $
            2.08000e-11,4.26000e-11,4.58000e-11,7.40000e-11]
; and deactivation by O2
chemq.q_81=[2.30000e-12,4.10000e-12,3.70000e-12,6.26000e-12,$
            4.90000e-12,5.19000e-12,4.10000e-12,2.60000e-12]

p_ave=[$
      [1.50979,      2.83820,      104.465,      16.9816],$
      [1.76115,      3.45162,      103.501,      17.1824],$
      [2.47268,      4.82708,      103.984,      16.9340],$
      [3.63904,      6.93513,      104.774,      16.5274],$
      [5.29935,      9.62446,      105.770,      16.0405],$
      [7.42349,      13.8895,      104.869,      16.4297],$
      [9.79957,      18.1602,      104.757,      16.4766],$
      [12.3029,      22.0481,      105.297,      16.2329] ]


chemq.q_92=5.e-17;N2p + N2    will never be important 
rateq.q_92 = chemq.q_92 * zminor.n_2p * zmaj.n2den

chemq.q_93=2.9e-11;N2p+NO
rateq.q_93 = chemq.q_93 * zminor.n_2p * zminor.no

chemq.q_96=9.5e-9;N2P + e -> N(4s,2d) + O (3P,1D);  Herron, ignoring t-dep
rateq.q_96 = chemq.q_96 * zminor.n_2p * zion.elec

;
; Radiative Recommbination
; ------------------------

; N(2D) -> N(4S) + hv   
chemr.r_5=fltarr(model.nlev) + 1.279e-5
rater.r_5 = chemr.r_5

chemr.r_94=fltarr(model.nlev) + 5.308e-3
rater.r_94 = chemr.r_94

chemr.r_95=fltarr(model.nlev) + 8.054e-2
rater.r_95 = chemr.r_95

;
; Branching Ratios
; ----------------
;

; ionization of O  ->  O+(2D) and O+(4S)   
F_O=0.66    

F_1S=.5
F_1D=.275

; it doesn't seem like this is labeled correctly....
; NO+ + e          ->  N(2D) and N(4S) and O   
F_2=.95;.85;5;5;from hellberg 2003

;this is the n2d yield from photodissociation of N2
;;THIS IS STOP GAP. CHANGE SOON;;
F_3S = 0.5; fltarr(model.nlev)
F_3D = 0.275
;F_3D = 0.0;fltarr(model.nlev)


; N2+ + O          ->  N(2D,4S) and NO+     R17
F_4=.90             
;    ->  O+ + N2 (1-F_4)
; Although it is clear that the N + NO+ channel is dominant, division among N(2D,4S)
; is less clear. Suggestive arguments for N2D as sole channel are given by dothe
; et al, JGR, 1996 and Mertens et al, GRl, 2008.  For N4S given by Scott et al, JPC,
; 1999 and measurements previous to Scott. Partitioning among N(2D)/N(4S) given by F_14


; N+ + O2 -> O2+ + N(2D,4S) 
F_5=.50 
;THis is the sum of the branching ratios for the 2 channels O2+ + N(2D) and O2+ + N(4S).  Division into (4S) and (2D) 
;effected by F_13.  Very convincing argument by Midey et al, JPCA, 2006 and  Dotan et al, Int.J.MassSpec.IonProc., 1997. 
;that N(2D) is the main atom product for this channel.
; See comment below F_6 
     
     
; N+ + O2          ->  NO+ + O
F_6=.42
;  Midey et al, JPCA, 2006.  Dotan et al, Int.J.MassSpec.IonProc., 1997
;Branching ratios F_5 and F_6 are new and update
;A third branch leading to products: O+ + NO
;has a branching of 1-F_5-F_6 = 0.06.  This change
;also meant several changes in the subroutines
;that calculate the production and losses, FGE 6-22-91

    
;  N+ + O2 -> O+ + NO    
F_7=.08
;  Midey et al, JPCA, 2006

;  N+ + NO -> NO+ + N(4S)
F_8=.91
;  Midey et al,JCP, 2004

;  N+ + NO -> N2+ + O 
F_9=.07
; Midey et al,JCP, 2004

; N+ + NO -> O+ + N2
F_10=.02
; Midey et al,JCP, 2004

;  N2A(v=3-6)+O->NO + N(2D). 
F_11=1.
; 1-F_11 is the BR for N(4S)

; N2+ +e->N(2D) 
F_12=.71
;from Petersen (1998)
;F_12p=.05;N2+ +e->N(2P) from Petersen (1998)
    
; N+ + O2 -> O2+ +N(4S,2D)    
F_13=1.     
;.5 ; F_5*F_13 gives the overall branching ratio of 
;N2D for R_25 while F_5*(1-F_13)  gives the same for N4S. F_13 is thus the percent of N(2D) in the O2+ + N
;channel. Very strong case by Midey et al, JPCA, 2006 and Dotan et al, Int.J.MassSpec.IonProc., 1997. that N(2D) is
; the main atom product for this channel, in which case F_13=1. 

;  N2+ + O -> N(2D,4S) + NO+   
F_14=.5;.5 
;As in F_13, F_14 is the percent of the N + NO+ channel going to the n4s channel
; Thus F_4 gives the actual N(2D)+NO+ branching raito, while
;(1 -F_14)*F_4) gives the same for N(4S).  Oran JGR (1975),Rusch (1991) and recent researchers
;(Mertens, GRL 2008) argue that N2D is the main channel

; N2D + O2-> NO + O(1D)
;         -> N(4S) + O2 
F_15=1.
; THere is a nonzero probability of quenching of the  N2D by the O2, rather than
; reaction. G_15 gives the total rate of removal of N2D by O2.  F_15 gives the amount
; that goes into the reactive channel.

;N2+e->n(2d,2p)+n(4s)
F_16=.26/.48
;  excited N br given by f1.  DIvision among the two states given by f_16.  i.e.  f_16=1
;means all the excited nitrogen is 2d. From Table of Zipf 1980.  After subtracting
;.01 from both N(4S) and N(2D) yields due to dissociative ionization, (included separately) the result is above. 


;N2p+O->N(2d,4s)+O
F_17=1.
;this gives the yield of n2d in the quenching of n(2p)
F_170=0.5 ; fraction of g_20 that leads to NO+

;n2p yield upon n2+ diss.rec.                        
F_18=.05
;from Petersen (1998)

;n2p yield from photodissociation
F_19=.05;
;.05; Fox 2007; n2p yield upon photodissociation (F_3 gives n2d yield)

;n2p yield for auroral dissociation (EPP)     
F_20=.05

chemf = {   F_O  : F_O ,$
            F_1s : F_1s,$
            F_1d : F_1d,$
            F_2  : F_2 ,$
            F_3s : F_3s,$
            F_3d : F_3d,$
            F_4  : F_4 ,$
            F_5  : F_5 ,$
            F_6  : F_6 ,$
            F_7  : F_7 ,$
            F_8  : F_8 ,$
            F_9  : F_9 ,$
            F_10 : F_10,$
            F_11 : F_11,$
            F_12 : F_12,$
            F_13 : F_13,$
            F_14 : F_14,$
            F_15 : F_15,$
            F_16 : F_16,$
            F_17 : F_17,$
            F_170: F_170,$
            F_18 : F_18,$
            F_19 : F_19,$
            F_20 : F_20 $
             }
             
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
