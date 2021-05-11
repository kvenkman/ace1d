;***************************************************************************************
;+
; NAME: mmddyyyyhhmm_d.pro
;
; LOCATION IN LIBRARY:
;
; TYPE:
;
; PURPOSE: Convert a date and time into a day number (1-366) and
;                  a fractional day remainder from the time
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUT ARGUMENTS: Variable Name         Type                    Description
;     MM = integer scalar or vector (1-12)
;     DD = integer scalar or vector (1-31)
;     YY = integer scalar or vector (00-99)
;          optional 4-digit value may be used
;     HH = integer scalar or vector (0-23), optional
;     MIN = integer scalar or vector (0-60), optional
;     SS = float scalar or vector (0-60), optional
;
; OUTPUTS: Variable Name         Type                    Description
;     floating point scalar or vector of day numbers
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; COMMON BLOCKS:
;
; CALLED ROUTINES:
;
; SIDE EFFECTS:
;
; RESTRICTIONS: This IDL routine is to be distributed only by the Naval
;				Research Laboratory, Code 7607.
;
; NOTES:
;
; MODIFICATION HISTORY:
;		Version 1.0: Scott Alan Budzien, NRL Code 7607, 10/4/97 Corrected for year 2000 bug
;-
;***************************************************************************************
function mmddyyyyhhmm_d,mm,dd,yy,hh,min,ss
;
; Parse input arguments
;
ny = n_elements(yy)
nm = n_elements(mm)
nd = n_elements(dd)
nh = n_elements(hh)
nmi = n_elements(min)
ns = n_elements(ss)
nmax = max([ny,nm,nd,nh,nmi,ns],listpnt)
list = ['YY','MM','DD','HH','MIN','SS']
format = "(a3,' must be either scalar or a vector of ',i0,' elements, e.g. "+$
         list(listpnt)+".')"
if (ny ne nmax) and (ny ne 1) then $
        message,string('YY',nmax,format=format)
if (nm ne nmax) and (nm ne 1) then $
        message,string('MM',nmax,format=format)
if (nd ne nmax) and (nd ne 1) then $
        message,string('DD',nmax,format=format)
if (nh ne 0) and (nh ne nmax) and (nh ne 1) then $
        message,string('HH',nmax,format=format)
if (nmi ne 0) and (nmi ne nmax) and (nmi ne 1) then $
        message,string('MIN',nmax,format=format)
if (ns ne 0) and (ns ne nmax) and (ns ne 1) then $
        message,string('SS',nmax,format=format)
;
; Create cumulative days previous to each month array
;  handle leap years
;
dpm = [-999,0,31,59,90,120,151,181,212,243,273,304,334]
;
; Convert possible scalar month, year values into an array
; Make sure it is compatible with the scalar or vector year value
;
if nm lt ny then mma = replicate(mm(0),ny) else mma = [fix(mm)]
if ny lt nm then yyyy = replicate(yy(0),nm) else yyyy = [fix(yy)]
;
; Perform leap year logic
;
mod4 = yyyy - yyyy/4*4
mod100 = yyyy - yyyy/100*100
mod400 = yyyy - yyyy/400*400
intercalary_day = (((mod4 eq 0) and (mod100 ne 0)) or (mod400 eq 0)) and (mma ge 3)
;
;
;
if nh eq 0  then hh = 0.
if nmi eq 0 then min = 0.
if ns eq 0  then ss = 0.
day = [double(dpm(mma)+dd)+hh/24.+min/1440.+ss/86400.] + intercalary_day
return,day
end
