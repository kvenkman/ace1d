PRO plot_diagnostics_helper, model_output, f107=f107, create_ps=create_ps
   
   ; load SNOE color table
   ; snoect, /bw
   a_env
   
   ; Model pressure coordinates
   zp = findgen(57)/4 - 7.
   
   !p.multi = [0, 2, 2, 0]
   if(keyword_set(create_ps)) then begin
      old_device = !d.name
      set_plot, 'ps'
      charthick = charthick
      ythick = 4
      xthick = xthick
      charsize = 1
      thick = 4
      
   endif else begin      
      window, /free, xsize = 800*(12/7.), ysize = 800
      charthick = 2
      thick = 2
      xthick = 2
      ythick = 2
      charsize=1.75
   endelse
   
   ; plot major species
   plot, alog10(model_output.zmaj.n2den), zp, yr = [-7, 7], title = 'Major Species', $
   background = 255, color = 0, charthick = charthick, xthick = xthick, $ 
   ythick = ythick, thick = thick, charsize = charsize, $
   xtitle = 'Log!D10!N (densities, cm!E-3!N)', ytitle = 'Z', /ys
   
   oplot, alog10(model_output.zmaj.o2den), zp, color = 3, thick = thick
   oplot, alog10(model_output.zmaj.oden), zp, color = 11, thick = thick

   xyouts, 12, 5.25, 'O', color = 11, charsize = charsize, charthick = charthick
   xyouts, 12, 4, 'O!D2!N', color = 3, charsize = charsize, charthick = charthick
   xyouts, 12, 2.75, 'N!D2!N', color = 0, charsize = charsize, charthick = charthick
   
   xyouts, 12, 7.5, 'F10.7:'+f107, charthick = charthick
   
   ; plot minor species
   plot, alog10(model_output.zminor.n4s), zp, yr = [-7, 7], xr = [0, 11], $
   title = 'Minor Species', background = 255, color = 0, charthick = charthick, $ 
   xthick = xthick, ythick = ythick, thick = thick, charsize = charsize, $
   xtitle = 'Log!D10!N (densities, cm!E-3!N)', ytitle = 'Z', /ys
   
   oplot, alog10(model_output.zminor.n4s_pce), zp, color = 0, thick = thick, linestyle = 2
   oplot, alog10(model_output.zminor.no_pce), zp, thick = thick, linestyle = 2, color = 3
   oplot, alog10(model_output.zminor.no), zp, color = 3, thick = thick
   oplot, alog10(model_output.zminor.n2d), zp, color = 11, thick = thick
   oplot, alog10(model_output.zminor.n2p), zp, color = 10, thick = thick
   oplot, alog10(model_output.zminor.o1d), zp, color = 6, thick = thick
   oplot, alog10(total(model_output.zminor.n2a, 1)), zp, color = 2, thick = thick
   
   xyouts, 9, 5.25, 'N(!E4!NS)', charsize = charsize,charthick = charthick, color = 255
   xyouts, 9, 4, 'N(!E2!ND)', color = 11, charsize = charsize,charthick = charthick
   xyouts, 9, 2.75, 'N(!E2!NP)', color = 10, charsize = charsize,charthick = charthick
   xyouts, 9, 1.5, 'NO', color = 3, charsize = charsize,charthick = charthick
   xyouts, 9, 0.25, 'O(!E1!ND)', color = 6, charsize = charsize,charthick = charthick
   xyouts, 9, -1.25, 'N!D2!N(A)', color = 2, charsize = charsize,charthick = charthick
   
   xyouts, 1950, 7.5, 'F10.7:'+f107, charthick = charthick, color = 255, charsize = 0.7
   
   ; plot temperatures
   plot, model_output.zmaj.tn, zp, xr = [100, 2200], background = 255, color = 0, $
   charthick = charthick, xthick = xthick, ythick = ythick, thick = thick, $
   charsize = charsize, xtitle = 'Temperature (K)', ytitle = 'Z', /ys

   oplot, model_output.zion.ti, zp, color = 3, thick = thick
   oplot, model_output.zion.te, zp, color = 11, thick = thick
   xyouts, 350, 5.25, 'Te', color = 11, charsize = charsize, charthick = charthick
   xyouts, 350, 4, 'Ti', color = 3, charsize = charsize, charthick = charthick
   xyouts, 350, 2.75, 'Tn', color = 0, charsize = charsize, charthick = charthick

   ; plot ionosphere
   plot, alog10(model_output.zion.e), zp, charthick = charthick, xthick = xthick, $
   xr = [1, 8], ythick = ythick, thick = thick, charsize = charsize, $
   xtitle = 'Log!D10!N (densities, cm!E-3!N)', ytitle = 'Z', $
   background = 255, color = 0, /ys
   ;title = 'Ion densities'
   oplot, alog10(model_output.zion.o_p), zp, color = 3, thick = thick
   oplot, alog10(model_output.zion.o2d_p), zp, color = 3, thick = thick, psym = 4
   oplot, alog10(model_output.zion.o2p_p), zp, color = 3, thick = thick, psym = 7
   oplot, alog10(model_output.zion.n_p), zp, color =  11, thick = thick
   oplot, alog10(model_output.zion.no_p), zp, color = 10, thick = thick
   oplot, alog10(model_output.zion.o2_p), zp, color = 6, thick = thick
   oplot, alog10(model_output.zion.n2_p), zp, color = 2, thick = thick

   xyouts, 7, 5.25, 'O!E+!N', color = 3, charsize = charsize, charthick = charthick
   xyouts, 7, 4, 'NO!E+!N', color = 10, charsize = charsize, charthick = charthick
   xyouts, 7, 2.75, 'O!D2!N!E+!N', color = 6, charsize = charsize, charthick = charthick
   xyouts, 7, 1.5, 'N!E+!N', color = 11, charsize = charsize, charthick = charthick
   xyouts, 7, .25, 'N!D2!N!E+!N', color = 2, charsize = charsize, charthick = charthick
   xyouts, 7, -1, 'e', charsize = charsize, charthick = charthick, color = 0

   plots, 6, -3.7, psym = 4, color = 3, thick = 4
   xyouts, 6.2, -4, 'O!E+!N(!E2!ND)', color = 3, charthick = charthick

   plots, 6, -4.95, psym = 7, color = 3, thick = 4
   xyouts, 6.2, -5.25,  'O!E+!N(!E2!NP)', color = 3, charthick = charthick

   xyouts, 6.5, 7.5, 'F10.7:'+f107, charthick = thick
  
   if(keyword_set(create_ps)) then begin
      device, /close
      set_plot, old_device
   endif

END

PRO plot_diagnostics, filename
   ; Given a ACE1D generated diagnostic save file, 
   ; plot model outputs
   
   if ~file_test(filename) then begin
      print, 'File: '+filename+" does not exist! Exiting."
      stop
   endif
   
   restore, filename
   
   ; Verify that the save file contains the correct output structure
   if ~keyword_set(model_output) then begin
      print, 'File does not contain ACE1D generated output structure. Exiting.'
      stop
   endif
   
   ; Verifying that the structure contains the expected fields
   names = tag_names(model_output)
   index1 = where(names EQ 'ACE1D_100')
   index2 = where(names EQ 'ACE1D_250')
   if (index1 EQ -1 || index2 EQ -1) then begin
      print, 'File was not generated by a diagnostic run. Exiting.'
      stop
   endif
   
   ; If the file name exists, contains the structure 'model_output', which 
   ; in turn contains the sub-structures ACE1D_100 and ACE1D_250, plot and 
   ; display model output fields
   
   plot_diagnostics_helper, model_output.ace1d_100, f107 = '100'
   ; plot_diagnostics_helper, model_output.ace1d_250, f107 = '250', /create_ps 
   
END
