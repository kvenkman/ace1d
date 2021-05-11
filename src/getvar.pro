function getvar,ncid,name,varinfo,help=help
if (keyword_set(help)) then begin
  print,format="(/,72('-'))"
  print,format="('getvar: read a variable from a currently open file')"
  print,format="(/,'Usage (called after a call to readfile)')"
  print,format="('  var = getvar(ncid,name,varinfo,[/help])')"
  print,format="(/,'Where:')"
  print,format="('  ncid is the current file id')"
  print,format="('  name is the name of the variable to read')"
  print,format="('  varinfo is a returned structure containing useful information about  the variable')"
  print,format="('  help is an optional keyword, to get this useage message.')"
  print,format="(/,'Example 1: Get this useage message')"
  print,format="('  h = getvar(/help)')"
  print,format="(/,'Example 2 (after calling readfile to get file structure): Read field TN')"
  print,'  TN = getvar(file.ncid,''TN'',vinfo)'
  print,' ' 
  print,'After the above call, execute ''help,vinfo,/struct'' to view the contents'
  print,'  of the returned structure.'
  print,format="(72('-'),/)"
  return,abort('')
endif
id = ncdf_varid(ncid,name)
if (id gt -1) then begin
; print,'Reading variable ',name,'...'
  ncdf_varget,ncid,id,f
  varinfo = ncdf_varinq(ncid,id)
  return,f
endif else begin
  print,'>>> WARNING: could not find variable with name "',name,'".'
  return,0
endelse
end