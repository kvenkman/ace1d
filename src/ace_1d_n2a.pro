PRO ace_1d_n2a, zminornow, zmajnow, zei, zminor, coeffmatrix, spindex, model
; Vibrational level specific N2(A) production rates are obtained by scaling
; the production of N2+ by photoelectrons, as described Yonker [2013]

; N2(A) losses are radiation in the VK band, quenching by and losses to O
; pages 38, 31, 83 (of document)
; g_v(z) = a + b*exp[-(z-c)^2/(2d^2)], where a - d are level specific fit parameters, z is the altitude

a = [1.50979, 1.76115, 2.47268, 3.63904, 5.29935, 7.42349, 9.79957, 12.3029]
b = [2.83820, 3.45162, 4.82708, 6.93513, 9.62446, 13.8895, 18.1602, 22.0481]
c = [104.465, 103.501, 103.984, 104.774, 105.770, 104.869, 104.757, 105.297]
d = [16.9816, 17.1824, 16.9340, 16.5274, 16.0405, 16.4297, 16.4766, 16.2329]

k_o =  [2.8, 3.4, 3.34, 3.48, 2.08, 4.26, 4.58, 7.4]*1e-11
k_o2 = [2.3, 4.1, 3.7, 6.26, 1,5.19, 4.1, 2.6]*1e-12
a_vk = [2.37, 2.65, 2.88, 3.10, 3.24, 3.32, 3.31, 3.31] ; Piper, 1993. Value for v=7 is assumed, but isn't important for NO/N(2D) calculations.
a_vk = 1./a_vk; Converting to loss frequencies

g = fltarr(8, model.nlev)
n2a_prod = fltarr(8, model.nlev) ; N2(A) production rate for v = 0-7
n2a_loss = fltarr(8, model.nlev) ; N2(A) loss frequency for v = 0-7
;n2a_prod_net = fltarr(model.nlev)
;n2a_loss_net = fltarr(model.nlev)
; coeffmatrix[spindex.n2a, spindex.o, *]
for i = 0, 7 do begin ; over vibrational levels
    g[i, *] = a[i] + b[i]*exp(-(zmajnow.zz-c[i])^2./(2.*d[i]^2.))
    n2a_prod[i, *] = zei.pei_n2*zmajnow.n2den/g[i, *]
    zminor.n2a[i, *] = n2a_prod[i, *]/(zmajnow.oden*k_o[i] + zmajnow.o2den*k_o2[i] + a_vk[i])
endfor

END
