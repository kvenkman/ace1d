FUNCTION ace_1d_trsolv, p, q, r, rhs


size1 = n_elements(p)

; Modified coefficient arrays
p_p = fltarr(size1)
q_p = fltarr(size1)
r_p = fltarr(size1)
rhs_p = fltarr(size1)

; Return value array
x = fltarr(size1)




r_p[0] = r[0]/q[0]
rhs_p[0] = rhs[0]/q[0]

FOR i = 1, size1 - 1 DO BEGIN
    r_p[i] = r[i]/(q[i] - p[i]*r_p[i-1])
    rhs_p[i] = (rhs[i] - p[i] * rhs_p[i-1])/(q[i] - p[i]*r_p[i-1])
ENDFOR

x[-1] = rhs_p[-1]
FOR i = size1 - 2, 0, -1. DO BEGIN
x[i] = rhs_p[i] - r_p[i]*x[i+1]
ENDFOR

return, x
END
