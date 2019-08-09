
function ace_1d_pcecalc,concin,prod,loss,dt

; this is a function to do a simple PCE calculation based on production and loss
; normally this is simply P/L, but if L=0, we want to avoid a divide by zero
; and just let the concentration be increased by P*dt

;ind=where(loss le 0.,n_ind);

;if n_ind eq 0.1 then begin
;
;  conc = prod / loss
;  
;endif else begin

  ;conc = conc + prod*dt 
  conc = (concin + prod*dt) / (1.0 + loss*dt)
  ;conc = concin + (prod*dt - loss*concin*dt)>0.0
  
;  if min(prod) lt 0. then stop,0,0
;  if max(prod) gt 1e8 then stop,0
;  if (1-min(finite(conc))) then stop,1
;  if (1-min(finite(prod))) then stop,2
;  if (1-min(finite(loss))) then stop,3
;endelse
  
return,conc
end

