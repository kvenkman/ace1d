
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
