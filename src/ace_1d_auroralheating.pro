PRO ace_1d_auroralheating, zmajnow, heatterms, model, pconst

; We assume an auroral energy input of 1keV particles with a flux
; of 0.025 ergs cm^-2 s^-1

flux = 1.56250e7 ; cm^-2 s^-1
epsilon = 1e2*pconst.q_e*1e7 ; eV -> J -> ergs

h = pconst.boltz*zmajnow.tn/(zmajnow.barm*pconst.grav/pconst.avo)
;rho = zmajnow.barm*model.p0*exp(-model.zp)/(pconst.boltz*zmajnow.tn*pconst.avo)
; We don't divide by rho as we want the heating rate in units of ergs cm^-3 s^-1
q_aurora = epsilon*flux/(3*h) ; K.D. Cole, 1975

heatterms.q_aurora = q_aurora ; We do not include this process
heatterms.q_total = heatterms.q_total + 0.*heatterms.q_aurora

; Smithtro [2005] does not include auroral electron heating 
; Roble [1987] states that solar min MSIS temperatures can be retrieved by only using an 
; auroral energy input of 13 GW, or by using 9.2 GW auroral energy and 70 GW Joule heating


END
