pro env_setup
   ;;; DISPLAY
   device, decomposed = 0, retain = 2

   ;;; set up the system variables
   !edit_input = 1000

   ct = { $
   black: 0, $
   violet: 1, $
   purple: 2, $
   navy: 3, $
   blue: 4, $
   aqua: 5, $
   green: 6, $
   lime: 7, $
   yellow: 8, $
   amber: 9, $
   orange: 10, $
   red: 11, $
   darkgrey: 100, $
   grey: 150, $
   lightgrey: 200, $
   white: 225}

   DEFSYSV, '!ct', ct, 1 ; 1 = read_only

   ;;; Set up the standard plotting
   !p.background=!ct.white
   !p.color=!ct.black
   !p.charsize=2
   !p.thick=3
   !p.charthick=3
   !x.thick=3
   !y.thick=3
   
   ; Load custom color table
   ; First 11 colors run through the spectrum, the rest are all grey
   ; 0 = black
   ; 1 = violet
   ; 2 = purple
   ; 3 = navy
   ; 4 = blue
   ; 5 = aqua
   ; 6 = green
   ; 7 = lime
   ; 8 = yellow
   ; 9 = amber
   ; 10 = orange
   ; 11 red
   ; 12 onward are all grey scale (dark to light)
   ; 255 white

   restore, "c11.tbl"
   tvlct,c1,c2,c3
   
end
