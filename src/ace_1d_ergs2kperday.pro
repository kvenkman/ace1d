function ace_1d_ergs2kperday,pconst,zmaj, hr

; from a given rate of reaction and exothermicity of that same reaction, 
; calculate a hearting rate in K per day

zmajn=(zmaj.oden+zmaj.o2den+zmaj.n2den) ; total number density of molecules
molecmass=zmaj.barm /pconst.avo *1e-3   ; mean mass of molecules, kg per molecule
heat = hr *1e6 / 1e7  /(zmajn*1e6*molecmass) / 1005. * (24.*60.*60.) ; Kelvin's per day, Cp is J K^-1 Kg^-1

return,heat
end
