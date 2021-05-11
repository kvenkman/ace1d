function new_model_atmos_par,model_name
;
; Create a generic parameter list that works for all atmospheric models
; --Additional parameters may be appended by the specific atmospheric
;   model code
;
time = met2ut(0)
if n_params(0) eq 0 then model_name = 'model'
par = create_struct(name=model_name+'_atm_par', 'model', model_name,$
                    time(0),$
                    'lat', 0., 'lon', 0., 'solar_local_time', 0., $
					'f10', 0., 'f10a', 0.)
return,par
end
