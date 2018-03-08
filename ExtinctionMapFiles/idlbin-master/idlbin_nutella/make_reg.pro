pro make_reg

restore, '/Users/jkrick/idlbin/object.sav'

headeracs = headfits('/Users/jkrick/hst/raw/wholeacs.fits')

a = where(object.acsm20 lt 0)
print, n_elements(a)
adxy, headeracs, object[a].ra, object[a].dec, x, y


openw, outlun, '/Users/jkrick/hst/raw/m20.reg', /get_lun
for i = 0l, n_elements(a) - 1 do begin
   printf, outlun, 'circle(', x(i),y(i), 20,')'
endfor

close, outlun
free_lun, outlun
end
