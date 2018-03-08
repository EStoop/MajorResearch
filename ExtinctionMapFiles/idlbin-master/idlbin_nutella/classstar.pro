pro classstar

close,/all
restore, '/Users/jkrick/idlbin/object.sav'        ;object
fits_read, '/Users/jkrick/hst/b_drz.fits', gdata, gheader
openw, outlun, '/Users/jkrick/hst/starlist', /get_lun
plothist, object.acsclassstar, xhist, yhist, /noprint, bin=0.01

for i = 0l, n_elements(object.acsclassstar) -1 do begin
   if object[i].acsclassstar gt 0.7 then begin
      adxy, gheader, object[i].ra, object[i].dec, x,y
      printf,outlun, x, y
   endif
endfor

close, outlun
free_lun, outlun


end
