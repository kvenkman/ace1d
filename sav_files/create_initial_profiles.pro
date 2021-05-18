PRO create_initial_profiles

   filenames = ["msis_080" , "msis_150", "msis_220"]
   
   for i = 0, n_elements(filenames) - 1 do begin
      ; we pre-allocate arrays. MSIS was run to generate altitude profiles between
      ; 70 and 1000 km with 5 km bins
      o = fltarr(187)
      n2 = fltarr(187)
      o2 = fltarr(187)
      tn = fltarr(187)
      
      OPENR, lun, filenames[i], /GET_LUN
      skip_lun, lun, 32, /lines ; Skip the header

      line_count = 0
      
      while not eof(lun) do begin
         readf, lun, alt, one, two, three, four
         
         o[line_count] = one
         n2[line_count] = two
         o2[line_count] = three 
         tn[line_count] = four
         
         line_count = line_count + 1
      endwhile
      print, "done reading ", filenames[i]
      free_lun, lun
      fit_to_model_grid, o, n2, o2, tn, filenames[i]
   endfor
END

PRO fit_to_model_grid, o, n2, o2, tn, filename
   ; Fit MSIS data to pressure grid which will be used by ACE1D
   
   ; Constants
   boltz = 1.38066e-16; Boltzmann's constant
   n_a = 6.023e23
   
   ; ACE1D model reference pressure and vertical coordinate grid
   model_p0 = 5e-4
   zp = findgen(57)/4 - 7.
   
   ; MSIS calculated number density, mass density, and pressure 
   n = o + n2 + o2  
   rho = (16.*o + 28.*n2 + 32.*o2)/n_a
   p = n*boltz*tn

   ; MSIS pressures expressed in the ACE1D coordinate system
   msis_p = -alog(p/model_p0)
   
   ; Interpolating number density to model grid   
   o_den  = interpol(o, msis_p, zp)
   o2_den = interpol(o2, msis_p, zp)
   n2_den = interpol(n2, msis_p, zp)
   tn = interpol(tn, msis_p, zp)

   ; Mass density (rho) interpolated to model grid
   rho_interp = interpol(rho, msis_p, zp)
   
   ; mass mixing ratios
   o_mmr = ((16.*o_den)/n_a)/rho_interp
   o2_mmr = ((32.*o2_den)/n_a)/rho_interp
   n2_mmr = ((28.*n2_den)/n_a)/rho_interp
   
   file_suffix = strmid(filename, 4, 4)
   
   save, o_den, o_mmr, o2_den, o2_mmr, n2_den, n2_mmr, tn, zp, filename = 'msis'+file_suffix+'.sav'
END
