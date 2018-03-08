pro rband
  close,/all
  restore, '/Users/jkrick/idlbin/object.sav'
  openw, outlun, '/Users/jkrick/hst/raw/rband.cat', /get_lun
  openw, outlun2, '/Users/jkrick/hst/raw/junk', /get_lun
;  fits_read, '/Users/jkrick/hst/raw/j9aj25010_drz.fits', data, header
  for i = 0l, n_elements(object.rmaga) - 1 do begin
     if object[i].rmaga gt 17 and object[i].rmaga le 24 and object[i].rfwhm lt 7.0 and object[i].rellip lt 0.2 then begin
;        adxy, header, object[i].ra, object[i].dec, x, y
        printf, outlun, object[i].ra, object[i].dec, object[i].rmaga
;        printf, outlun2, x, y
     endif
  endfor
  close, /all
 
end
