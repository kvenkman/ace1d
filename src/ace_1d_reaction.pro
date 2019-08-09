
pro ace_1d_reaction,index_r1,index_r2,index_p1,index_p2,rate,exo,chemmat,heatmat,spname,reaction,yield_p1=yield_p1

; 
; routine to adjust chemistry and heat matrices for a givenr eaction
; note that quenching reactions have a reactant and product that is the same - need to make sure this is handled right in P and L calculations
; 
; 
; we need to make sure we populate arrays with higher index first
if index_r1 gt index_r2 and (index_r2 ne -1) then begin
ir1 = index_r1
ir2 = index_r2
endif
if index_r2 gt index_r1 then begin
ir1 = index_r2
ir2 = index_r1
endif

if (not keyword_set(yield_p1)) then yield_p1 = 1.0

; product one must always be supplied

; reaction rate placed into production / loss matrix
chemmat[ir1,ir2,index_p1] = chemmat[ir1,ir2,index_p1] + rate * yield_p1
chemmat[ir1,ir2,index_p2] = chemmat[ir1,ir2,index_p2] + rate

; chemical heating (do only once!)
heatmat[ir1,ir2,index_p1] = heatmat[ir1,ir2,index_p1] + ace_1d_rate2heat(pconst,zmaj,rate,exo)

; name the reaction for plotting later
; watch out for 1body reactions or single products
if      spname[ir2] eq '1Body' then s2 = '' else s2 = ' + ' + spname[ir2]
if spname[index_p2] eq '1Body' then s3 = '' else s3 = ' + ' + spname[index_p2]

reaction[ir1,ir2,index_p1] = spname[ir1] + s2 + ' --> ' + spname[index_p1] + s3

return
end



