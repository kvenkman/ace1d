
pro ace_read_ftable,lun,x0,table,label

n=43

label='a string'
s=label
readf,lun,label

x0=fltarr(n)
table=fltarr(n,9)
dum=fltarr(9)

for i=0,n-1 do begin

readf,lun,a,dum
x0[i]=a
table[i,*]=dum

endfor

return
end