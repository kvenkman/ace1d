;-------------------------------------------------------------
;+
; NAME:
;       FJULIAN
; PURPOSE:
;       From Year, Month, and Day compute Julian Day number.
;       This routine properly handles fractional days, unlike the
;       built-in IDL or APL routines.
; CATEGORY:
; CALLING SEQUENCE:
;       jd = FJULIAN(y,m,d)
; INPUTS:
;       y = Year (like 1987).                    in
;       m = month (like 7 for July).             in
;       d = month day (like 23).                 in
; KEYWORD PARAMETERS:
;       /HELP    print out parameter list
;       /MJD     modified Julian date (days since 2400000.5, or Nov 17, 1858 0 UT)
; OUTPUTS:
;       jd = Julian Day number (like 2447000).   out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner,  23 June, 1985 --- converted from FORTRAN.
;       Johns Hopkins University Applied Physics Laboratory.
;       RES 18 Sep, 1989 --- converted to SUN
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------

	function fjulian, iy, im, id, help=hlp, mjd=mjd

;
; Parse function call
;
	if (n_params(0) LT 3) or keyword_set(hlp) then begin
	  print,' From Year, Month, and Day compute floating point Julian Day number.'
	  print,' jd = fjulian(y,m,d)'
	  print,'   y = Year (like 1987).       '
	  print,'   m = month (like 7 for July).'
	  print,'   d = month day (like 23).    '
      print,'   /MJD = optional modified Julian date since 24400000.5'
	  return, -1
	endif
;
; Parse parameters
;
    ny = n_elements(iy)
    nm = n_elements(im)
    nd = n_elements(id)
    n = max([ny,nm,nd])
    if (ny ne 1 and ny ne n) or $
       (nm ne 1 and nm ne n) or $
       (nd ne 1 and nd ne n) then begin
        message,"Mismatched parameter array sizes.",/info
        return,-1
    endif
;
; convert parameters adn make sure they are either scalar or vector
;
	if ny eq 1 then y = long(iy(0)) else y = long(iy)
	if nm eq 1 then m = long(im(0)) else m = long(im)
	if nd eq 1 then d = double(id(0)) else d = double(id)
;
; calculation
	jd = 367*y-7*(y+(m+9)/12)/4-3*((y+(m-9)/7)/100+1)/4 $
             +275*m/9+d+1721028.5d0
;
; Modified Julian date option
;
    if keyword_set(mjd) then jd = jd - 2400000.5d0
	return, jd

	end
