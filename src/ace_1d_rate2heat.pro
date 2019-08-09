function ace_1d_rate2heat,pconst,zmaj, rate,exo

; from a given rate of reaction and exothermicity of that same reaction, 
; calculate a hearting rate in K per day

zmajn=(zmaj.oden+zmaj.o2den+zmaj.n2den) ; total number density of molecules
molecmass=zmaj.barm /pconst.avo *1e-3   ; mean mass of molecules, kg per molecule
heat = rate * exo * 1e6 * !a.conv.j_per_ev /(zmajn*1e6*molecmass) / pconst.air_cp * (24.*60.*60.) ; Kelvin's per day, Cp is J K^-1 Kg^-1

return,heat
end
