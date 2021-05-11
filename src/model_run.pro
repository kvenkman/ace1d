PRO model_run, f107, outputs

   ; this parameter can be changed to mask out certain 
   ; solar wavelength bins
   mask = fltarr(37) + 1. 
    
   ; Set the indices of the bins to be zeroed out here
   mask_idx = [] 
   mask[mask_idx] = 0
	
   inputs = { $ 
      run_year : 1999 , $
      start_day: 80   , $
      ndays    : 2.*15, $
      runlat   : 0    , $
      timestep : 60.  , $
      save_res : 2*24l, $
      f107d    : f107 , $
      f107a    : f107 , $
      mask     : mask   $
   }

   ace_1d_mainprogram, inputs, outputs

   print, "Done with F10.7: ", f107
   return	
END

   ;; The code below calls the model_run procedure, which in turn runs the ACE1D code
   output_filename = 'solar_min_max_run.sav'

   f107 = [100, 250]

   save_exe_str = "model_output = {"
   for i = 0, n_elements(f107) - 1 do begin
      exe_str = "model_run, " + strtrim(f107[i], 2) + ", ace1d_"+strtrim(f107[i], 2)
      value = execute(exe_str)
      save_exe_str = save_exe_str + "ace1d_"+strtrim(f107[i], 2)+":"+"ace1d_"+strtrim(f107[i], 2) + ","
   endfor

   save_exe_str = strmid(save_exe_str, 0, strlen(save_exe_str)-1) + "}"

   value =	execute(save_exe_str)

   save, model_output, filename = output_filename

END

