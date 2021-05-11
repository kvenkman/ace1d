; restore, 'model_output_tn_70gw_n2drate_test.sav'
; restore, 'model_output_tn_70gw.sav'
; restore, 'model_output_tn_70gw_eddy_test.sav'

restore, 'model_output_tn_70gw.sav'
tn1 = t_exo
restore, 'model_output_tn_70gw_zero_nov.sav'
tn2 = t_exo

diff = (tn2-tn1)/tn2
print, 100.*diff


END
