pro yd_mmddyyyyhhmmss,yd_in,mm,dd,yyyy,hr,mi,ss
;+
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; YD_MMDDYYYYHHMM:  Convert a Year/day into a date and time.
;
; written by: Scott Alan Budzien
; completed :
;
; revision history:
; 1)  Corrected for year 2000 bug 10/4/97
; 2)
; 3)
; 4)
;
; input arguments:
;     YD = YYDDD.fractional day or YYYYDDD.fractional day
;
; output:
;     scalars or vector of date/time parameters
;
; output keywords:
;     None.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;-
;
; Parse input arguments
;
nmax = n_elements(yd_in)
yyyy = long(yd_in)/1000
yd = yd_in - yyyy*1000
iday = fix(yd)
fday = yd-iday
;
; Create cumulative days previous to each month array
;  handle leap years
;
;      0  1	 2	3	4	5	6	7	8	9  10  11
dpm = [0,31,59,90,120,151,181,212,243,273,304,334]
dpmleap = dpm
dpmleap(2:*) = dpmleap(2:*) + 1
;
; Perform leap year logic
;
mod4 = yyyy - yyyy/4*4
mod100 = yyyy - yyyy/100*100
mod400 = yyyy - yyyy/400*400
intercalary_day = ((mod4 eq 0) and (mod100 ne 0)) or (mod400 eq 0)
;
; Convert YYYY, iday to DATE
;
if nmax eq 1 then begin
    if intercalary_day(0) then begin
        sel = max(where(iday gt dpmleap))
        mm = sel(0)+1
        dd = iday - dpmleap(sel)
    endif else begin
        sel = max(where(iday gt dpm))
        mm = sel(0)+1
        dd = iday - dpm(sel)
    endelse
endif else begin
    mm = intarr(nmax)
    dd = intarr(nmax)
    for i=0,nmax-1 do $
        if intercalary_day(i) then begin
            sel = max(where(iday(i) gt dpmleap))
            mm(i) = sel+1
            dd(i) = iday(i) - dpmleap(sel)
        endif else begin
            sel = max(where(iday(i) gt dpm))
            mm(i) = sel+1
            dd(i) = iday(i) - dpm(sel)
        endelse
endelse
;
; convert fractional day to time
;
hr = fix(fday*24)
mi = fix(fday*1440L) - hr*60
ss = float(fday*86400d0 - hr*3600d0 - mi*60d0)
return
end
