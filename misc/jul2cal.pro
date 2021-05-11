function jul2cal, Julian_day, Frac_Julian_day
; NAME:
;   JUL2CAL
;
; PURPOSE:
;   Return the calendar date and time given julian date.
;   This is the inverse of the function FJULIAN.  Note that
;   these custom-designed functions properly handle the 0.5
;   day offset between Julian and Solar days, whereas the
;   built-in IDL routines do not!
;
; CATEGORY:
;   Misc.
;
; CALLING SEQUENCE:
;   out = CALDAT(Julian,Frac_Julian)
;   See also: fjulian, the inverse of this function.
;
; INPUTS:
;   JULIAN contains the Julian Day Number (which begins at noon) of the
;   specified calendar date.  It can be a long integer or double float.
;   FRAC_JULIAN contains the fractional julian day (which begins at noon) of the
;   specified calendar date.  It should be a float or double.
;
; OUTPUTS:
;   Structure containing
;   (Trailing parameters may be omitted if not required.)
;   YR: Number of the desired year.
;   MO: Number of the desired month (1 = January, ..., 12 = December).
;   DY: Number of day of the month, including fractional extension.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   None.
;
; RESTRICTIONS:
;   Accuracy using IEEE double precision numbers is approximately
;   1/10000th of a second.
;
; MODIFICATION HISTORY:
;   Translated from "Numerical Recipies in C", by William H. Press,
;   Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
;   Cambridge University Press, 1988 (second printing).
;
;   DMS, July 1992.
;   DMS, April 1996, Added HOUR, MINUTE and SECOND keyword
;   SAB, Sep 2, 1998, Generalized to handle array input.
;-
;
;
; Parse parameters
;
if n_params(0) eq 2 then begin
    n1 = n_elements(julian_day)
    n2 = n_elements(frac_julian_day)
    if n1 eq 1 and n2 eq 1 then begin
        julian = long(julian_day(0))
        frac_julian = double(frac_julian_day(0))
        n = 1
    endif else if n1 eq 1 and n2 gt 1 then begin
        julian = replicate(long(julian_day(0)),n2)
        n = n2
    endif else if n1 gt 1 and n2 eq 1 then begin
        frac_julian = replicate(double(frac_julian_day(0)),n1)
        n = n1
    endif else if n1 eq n2 then begin
        julian = long(julian_day)
        frac_julian = double(frac_julian_day)
        n = n1
    endif else begin
        message,"Mismatched or invalid parameter array sizes.",/info
        return,-1
    endelse
endif else if n_params(0) eq 1 then begin
    n = n_elements(julian_day)
    julian = long(julian_day)
    frac_julian = double(julian_day-julian)
endif
;
IGREG = 2299161L    ;Beginning of Gregorian calendar

julian0 = long(Julian)
frac = Frac_Julian - 0.5
s = where(frac lt 0,count)
if count gt 0 then frac(s) = frac(s)+1.

s = where(frac lt 0.5,count)
if count gt 0 then julian0(s) = julian0(s) + 1

ja = julian0
s = where(julian0 ge igreg,count)
if count gt 0then begin
    jalpha = long(((julian0(s) - 1867216l) - 0.25d0) / 36524.25d0)
    ja(s) = julian0(s) + 1 + jalpha - long(0.25d0 * jalpha)
endif

jb = ja + 1524
jc = long(6680.0d0 + ((jb-2439870l)-122.1d0)/365.25d0)
jd = long(365 * jc + (0.25d0 * jc))
je = long((jb - jd) / 30.6001d0)

day = jb - jd - long(30.6001 * je)
month = je -1
year = jc - 4715

s = where(month gt 12,count)
if count gt 0 then month(s) = month(s) - 12

s = where(month gt 2,count)
if count gt 0 then year(s) = year(s) - 1

s = where(year le 0,count)
if count gt 0 then year(s) = year(s) - 1
;
; Create output structure
;
s0 = {yr:0,mo:0,dy:0d}
out = replicate(s0,n)
if n_elements(year) gt 1 then out.yr = reform(year) else out.yr = year
if n_elements(month) gt 1 then out.mo = reform(month) else out.mo = month
if n_elements(day+frac) gt 1 then out.dy = reform(day+frac) else out.dy = day+frac
return,out
end
