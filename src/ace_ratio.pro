
function ace_ratio,fnum,fden

; this procedure simply takes the ratio fnum / fden. If an element of den is 0.0 though, a value of zero is returned
; the goal is to avoid 1./0. = NaN issues

ratio=fnum-fnum                                                   ; make array of all zeroes
indne0=where(fden gt 1e-20,n_indne0)                                ; find all locations where a ratio is valid
if n_indne0 gt 0 then ratio[indne0] = fnum[indne0] / fden[indne0] ; calc ratio at valid places

return, ratio
end
