
; Running the ACE 1D for various levels of solar actvity:

pro runonce,sa,output

inputs = {run_year:1999, start_day:80, ndays:14, runlat:0, timestep:120., $
 save_res:2*24l,f107d:sa,f107a:sa} ; save resolution is given in minutes
 ace_1d_mainprogram,inputs,output
 
return
end

runonce,70.,ace1d_70
runonce,100.,ace1d_100
runonce,150.,ace1d_150
runonce,175.,ace1d_175
runonce,200.,ace1d_200
runonce,210.,ace1d_210
runonce,220.,ace1d_220
runonce,230.,ace1d_230
runonce,240.,ace1d_240
runonce,245.,ace1d_245
runonce,250.,ace1d_250
 
save,file='ace1d_solaractivity_halfo3p_zeronov.sav',ace1d_70,ace1d_100,ace1d_150,ace1d_175,$
ace1d_200,ace1d_210,ace1d_220,ace1d_230,ace1d_240,ace1d_245,ace1d_250

end

