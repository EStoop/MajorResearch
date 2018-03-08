pro massvscolor

restore, '/Users/jkrick/idlbin/object.sav'
ps_open, filename='/Users/jkrick/nep/massvscolor.ps',/portrait,/square,/color
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader


plothist, alog10(object.mass), xhist, yhist, bin=0.1, /noprint,/nan;,/noplot;,xrange=[0,480]

plot, object.gmaga, object.gfwhm, psym = 2, xrange=[5,30], yrange=[0,30], xtitle='gmag',ytitle = 'g fwhm'
plot, object.gmaga, object.gellip, psym = 2, xrange=[5,30], yrange=[0,1], xtitle='gmag',ytitle = 'g ellip'

;what are the objects with gfwhm lt 7?  stars and galaxies


for i = 0, n_elements(object.ra) -1 do begin
;   if object[i].gfwhm lt 7 and object[i].gfwhm gt 0 then begin
;      adxy, gheader, object[i].gxcenter, object[i].gycenter,x,y
;      print, x, y
;   endif
;   if object[i].gmaga lt 18 and object[i].gmaga gt 0 and object[i].gfwhm lt 20 then begin
;      adxy, gheader, object[i].gxcenter, object[i].gycenter,x,y
;      print, x, y
;   endif
   if object[i].mass gt 1E12 then begin
      adxy, gheader, object[i].gxcenter, object[i].gycenter,x,y
      print, x, y
   endif
endfor



u = fltarr(n_elements(object.umaga))
g = fltarr(n_elements(object.gmaga))
j = 0
k = 0



for i = 0, n_elements(object.ra) - 1 do begin
   if object[i].umaga gt 10 and object[i].umaga lt 30 then begin
      u[j] = object[i].umaga
      j = j + 1
   endif
   if object[i].gmaga gt 10 and object[i].gmaga lt 30 then begin
      g[k] = object[i].gmaga
      k = k + 1
   endif


endfor

u = u[0:j-1]
g = g[0:k-1]

plot, alog10(object.mass), u-g,psym=3,xrange=[6,15],yrange=[-3,3]
plot, alog10(object.mass), object.umaga-object.gmaga,psym=3,xrange=[6,15],yrange=[-3,3]

ps_close, /noprint,/noid
end
