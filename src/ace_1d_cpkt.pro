PRO ace_1d_cpkt, zmaj, mass, pconst, model

; This procedure defines specific heat capacity and thermal conductivity coefficients
; Ref Aeronomy Banks & Kockarts, part B. pg. 4, 12-14

      vmr_o2 = zmaj.barm*zmaj.o2/mass.o2
      vmr_o  = zmaj.barm*zmaj.o/mass.o
      vmr_n2 = 1. - vmr_o -vmr_o2 ;zmaj.barm*zmaj.n2/mass.n2

      zmaj.cp = pconst.gask*.5*(vmr_o2*7./32.+vmr_n2*7./28.+ $
                                              vmr_o*5./16.)
                                               
      zmaj.kt = (56.*(vmr_o2 + vmr_n2) + 75.9*vmr_o)*((zmaj.tn)^0.69)
      
END
