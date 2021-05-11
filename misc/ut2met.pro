;***************************************************************************************
;+
; NAME: ut2met.pro
;
; LOCATION IN LIBRARY:
;
; TYPE:
;
; PURPOSE:
; 	Converts UT values from missions including ARGOS, STS-39 (UVLIM)
; 	mission into a simple MET value in mission-defined "ticks" (usually
;       seconds).
;
; Special Note:  UVLIM MET is generally stated in deciseconds!
;
; CATEGORY:
;
; CALLING SEQUENCE: ut2met,yyyy,mm,dd,mission_name
;
; INPUT ARGUMENTS: Variable Name         Type                    Description
;                   yyyy                 long                    four-digit year
;                     mm                 int                     month number 1-12
;                     dd                 double                  date+fractional day (e.g. 1.000-31.9999999...)
;
; OUTPUTS: Variable Name         Type                    Description
;	return		         double			MET in "ticks"
;
; INPUT KEYWORDS: these are included for backwards compatability
;	argos		argos P91 mission, this is the default
;	sts39		sts-39 UVLIM mission
;	STRV1D		CERTO-PLUS mission
;	f16   		DMSP F16 mission (SSULI)
;	f17   		DMSP F17 mission (SSULI)
;       COSMIC		FORMOSAT3/COSMIC (TIP)
;       ANDERR		STS-116 ANDE Risk Reduction Deploy
;   f18          DMSP F17 mission (SSULI)
;       COSMIC   FORMOSAT3/COSMIC (TIP)
;       ANDERR   STS-116 ANDE Risk Reduction
;       ANDE     STS-127 ANDE
;       RAIDS    HTV-1 RAIDS to ISS
;
; OUTPUT KEYWORDS:
;
; COMMON BLOCKS:
;       none, but requires user-defined system variable
;
; CALLED ROUTINES:
;	fjulian
;	launchtimesinit
;
; SIDE EFFECTS:
;	may define system variable !launchtimes by running launchtimesinit if not defined
;
; RESTRICTIONS: This IDL routine is to be distributed only by the Naval
;				Research Laboratory, Code 7607.
;
; NOTES:
;
; MODIFICATION HISTORY:
;		Version 0.9: Scott Alan Budzien, NRL Code 7607, 07/07/99  met2ut
;		Version 1.1: Clyde Fortna, SFA inc. @NRL Code 7607, 07/17/2002 split off launchtimes structure
;-
;***************************************************************************************
function ut2met,yyyy,mm,dd,mission_name,argos=argos,sts39=sts39,STRV1D=strv1d,$
         F16=f16,F17=f17,F18=f18,cosmic=COSMIC,anderr=ANDERR,ande=ANDE,raids=RAIDS
;
; backward compatibility
;
CASE n_params() OF
    3: BEGIN
       ; set the mission string from keywords
       if      keyword_set(sts39) then mission_name='STS-39' $
       else if keyword_set(strv1d) then mission_name='STRV-1D' $
       else if keyword_set(f16) then mission_name='F16' $
       else if keyword_set(cosmic) then mission_name='COSMIC' $
       else if keyword_set(f17) then mission_name='F17' $
       else if keyword_set(f18) then mission_name='F18' $
       else if keyword_set(anderr) then mission_name='ANDERR' $
       else if keyword_set(ande) then mission_name='ANDE' $
       else if keyword_set(raids) then mission_name='RAIDS' $
       else mission_name='ARGOS' ; default
       END
    4: mission_name = strupcase(mission_name)
    else: begin
          message,"Requires year, month, date and optional mission name as parameters."
          end
ENDCASE
;
; Make sure system parameter has been defined
;
DEFSYSV, '!launchtimes', EXISTS = exists
IF NOT exists THEN launchtimesinit ; initialize system value if not already there
;
; match mission name to index in array
;
s=(where(!launchtimes.mission_name eq mission_name))(0)
if s lt 0 then $
	message, "mission_name must be one of:'+!launchtimes.mission_name
;
; Calculate launch time Julian date
;
launch_sec = !launchtimes(s).gpssecofweek
launch_week = !launchtimes(s).gpsweek
launch_jd = fjulian(!launchtimes(s).year,!launchtimes(s).month, $
                    !launchtimes(s).day + (launch_sec mod 86400)/86400d0 )
;
; Calculate parameter Julian date
;
jd = fjulian(yyyy,mm,dd)
;
; Don't forget to express the output in "ticks"
;
return,(jd-launch_jd)*!launchtimes(s).ticksperday
end
