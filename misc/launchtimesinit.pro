;***************************************************************************************
;+
; NAME: launchtimesinit.pro
;
; LOCATION IN LIBRARY:
;
; TYPE:
;
; PURPOSE:
;	Initialize the system variable which shall store the launch time UT structure
;	for each space mission for converting between UT and MET.
;	UT is universal time and MET is mission elapsed time in seconds with possible fraction.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUT ARGUMENTS: Variable Name         Type                    Description
;		dummy		n/a		not used
;
; OUTPUTS: Variable Name         Type                    Description
;		fills system variable !launchtimes with the mission launch time structure
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
;		Version 0.9: Scott Alan Budzien, NRL Code 7607, 07/07/99  met2ut
;		Version 1.1: Scott Alan Budzien, NRL Code 7607, 07/17/2002 split into launchtimes structure
;-
;***************************************************************************************
pro launchtimesinit,dummy
;
; Number of NRL launches:
n_launches = 10
;
lt0 = {launchtimestruct,mission_name:'',$
       year:0,$
       month:0,$
       day:0,$
       hour_ut:0,$
       minute:0,$
       second:0.,$
       gpsweek:0,$
       gpssecofweek:0d,$
       ticksperday:0D}
launchtimes = replicate(lt0,n_launches)
;
; mission index, mi
mi = 0    
;
; Define constants
;
launchtimes(mi).mission_name  = 'STS-39'
launchtimes(mi).year          = 1991
launchtimes(mi).month         = 4
launchtimes(mi).day           = 28
launchtimes(mi).hour_ut       = 11
launchtimes(mi).minute        = 33
launchtimes(mi).second        = 15.
launchtimes(mi).gpsweek       = 590
launchtimes(mi).gpssecofweek  = 41595L      ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 864000d0
;
launchtimes(mi).mission_name  = 'ARGOS'
launchtimes(mi).year          = 1999
launchtimes(mi).month         = 2
launchtimes(mi).day           = 23
launchtimes(mi).hour_ut       = 10
launchtimes(mi).minute        = 29
launchtimes(mi).second        = 55
launchtimes(mi).gpsweek       = 998
launchtimes(mi).gpssecofweek  = 210595L      ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0
;
launchtimes(mi).mission_name  = 'STRV-1D'
launchtimes(mi).year          = 2000
launchtimes(mi).month         = 11
launchtimes(mi).day           = 16
launchtimes(mi).hour_ut       = 1
launchtimes(mi).minute        = 07
launchtimes(mi).gpsweek       = 1088
launchtimes(mi).gpssecofweek  = 349620L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0
;
launchtimes(mi).mission_name  ='F16'
launchtimes(mi).year          = 2003
launchtimes(mi).month         = 10
launchtimes(mi).day           = 18
launchtimes(mi).hour_ut       = 16
launchtimes(mi).minute        = 17
launchtimes(mi).second        = 0
launchtimes(mi).gpsweek       = 1240
launchtimes(mi).gpssecofweek  = 577020L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0
;
launchtimes(mi).mission_name  = 'COSMIC'
launchtimes(mi).year          = 2006
launchtimes(mi).month         = 4
launchtimes(mi).day           = 15
launchtimes(mi).hour_ut       = 1
launchtimes(mi).minute        = 40
launchtimes(mi).second        = 0
launchtimes(mi).gpsweek       = 1370
launchtimes(mi).gpssecofweek  = 524400L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0
;

launchtimes(mi).mission_name  = 'F17'

launchtimes(mi).year          = 2006

launchtimes(mi).month         = 11

launchtimes(mi).day           = 4

launchtimes(mi).hour_ut       = 13

launchtimes(mi).minute        = 53

launchtimes(mi).second        = 02

launchtimes(mi).gpsweek       = 1399

launchtimes(mi).gpssecofweek  = 568382L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)

launchtimes(mi++).ticksperday = 86400d0

;

launchtimes(mi).mission_name  = 'ANDERR'

launchtimes(mi).year          = 2006

launchtimes(mi).month         = 12

launchtimes(mi).day           = 21

launchtimes(mi).hour_ut       = 18

launchtimes(mi).minute        = 23

launchtimes(mi).second        = 15

launchtimes(mi).gpsweek       = 1406

launchtimes(mi).gpssecofweek  = 411795L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)

launchtimes(mi++).ticksperday = 86400d0
;
launchtimes(mi).mission_name  = 'ANDE'
launchtimes(mi).year          = 2009
launchtimes(mi).month         = 7
launchtimes(mi).day           = 15
launchtimes(mi).hour_ut       = 22
launchtimes(mi).minute        = 3
launchtimes(mi).second        = 10
launchtimes(mi).gpsweek       = 1540
launchtimes(mi).gpssecofweek  = 338592L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0

;
launchtimes(mi).mission_name  = 'RAIDS'     ; 2009-09-11 02:01:46 JST (+9 GMT)
launchtimes(mi).year          = 2009
launchtimes(mi).month         = 9
launchtimes(mi).day           = 10
launchtimes(mi).hour_ut       = 17
launchtimes(mi).minute        = 1
launchtimes(mi).second        = 46
launchtimes(mi).gpsweek       = 1548
launchtimes(mi).gpssecofweek  = 406906L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0
;
launchtimes(mi).mission_name  = 'F18'
launchtimes(mi).year          = 2009
launchtimes(mi).month         = 10
launchtimes(mi).day           = 18
launchtimes(mi).hour_ut       = 16
launchtimes(mi).minute        = 12
launchtimes(mi).second        = 0
launchtimes(mi).gpsweek       = 1554
launchtimes(mi).gpssecofweek  = 58320L     ;this is NOT intrinsic GPS time!  (leap secs have been adjusted)
launchtimes(mi++).ticksperday = 86400d0
;
defsysv,'!launchtimes',launchtimes
;
return
end
