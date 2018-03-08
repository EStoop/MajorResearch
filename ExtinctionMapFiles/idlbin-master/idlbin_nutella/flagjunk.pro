pro flagjunk
close, /all

fits_read, '/Users/jkrick/hst/raw/wholeacs.fits', data, header
restore, '/Users/jkrick/idlbin/object.sav'
openw, outlun, '/Users/jkrick/hst/raw/trashpoints.reg',/get_lun

targ = where (object.acsmag gt 22 and object.acsmag lt 90.)
;targ = where (x lt 4300. and x gt 1400. and y gt 6700 and y lt 10000)
print, n_elements(targ)
c = 0
adxy, header, object[targ].acsra, object[targ].acsdec, x, y
aper, data, x, y, flux, err, sky, skyerr, 2, [25], [35,45], [1,1],/flux,/silent

for i = 0l, n_elements(targ) - 1 do begin
;   print, fluxaper(targ(i)), flux,  sky, skyerr,format='(F10.2, F10.2, F10.2, F10.2)'
   thresh = 6.0*1.5* sqrt(sky(i))
;   print, thresh
;   if thresh gt fluxaper(targ(i)) then printf, outlun, 'circle(', x(targ(i)),y(targ(i)), 50,')'
   if thresh gt flux(i) and object[targ(i)].rmaga gt 90 and   object[targ(i)].irac1mag gt 90  then begin ;or thresh gt 10^((25.937 - object[targ(i)].acsmag)/(-2.5)) then begin
;      print, sky(i), skyerr(i), flux(i), thresh, 10^((25.937 - object[targ(i)].acsmag)/(-2.5)), x(i), y(i), format='(F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2)'
      print, object[targ(i)].rmaga, object[targ(i)].irac1mag
      printf, outlun, 'circle(', x(i),y(i), 30,')'
      c = c + 1
   endif

endfor

print, c
close, outlun
free_lun, outlun
end


;  readcol,'/Users/jkrick/palomar/lfc/catalog/test.cat',num, x, y, fluxauto, magauto, magerrauto, fluxaper, magaper, magerraper, fwhm, isoarea, back, thresh, mu_thresh,format="A"
