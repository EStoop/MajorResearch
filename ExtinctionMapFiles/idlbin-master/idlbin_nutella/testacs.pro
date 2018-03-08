pro test_acs

close,/all
restore, '/Users/jkrick/idlbin/object.sav'

ps_open, filename='/Users/jkrick/hst/phot.ps',/portrait,/square,/color

;----------------------------------------------------------------
 readcol, '/Users/jkrick/hst/test.ASC', id, ra, dec, x,y,magiso, magisocor, magaper, magauto, magerrauto, magbest, kronrad, back, thresh, fluxmax,isoarea, theta, flag,fwhm,elong,ellip ,classstar, format="A"

print, ra(2), dec(2)

; create initial arrays

iarr = fltarr(n_elements(id))
magisoarr= fltarr(n_elements(id))
magisocorarr= fltarr(n_elements(id))
magaperarr= fltarr(n_elements(id))
magautoarr= fltarr(n_elements(id))
magbestarr= fltarr(n_elements(id))
fwhmarr =  fltarr(n_elements(id))
xarr = fltarr(n_elements(id))
yarr = fltarr(n_elements(id))
raarr = fltarr(n_elements(id))
decarr = fltarr(n_elements(id))

countarr = fltarr(n_elements(id))
deltaarr = fltarr(n_elements(id))

m=n_elements(id)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999

cind = 0
  
dist=irmatch
dist[*]=0
openw, outlun2, '/Users/jkrick/hst/phot.out', /get_lun
for q=long(0),m-1 do begin
   dist=sphdist( ra(q), dec(q),object.rxcenter,object.rycenter,/degrees)

  ;count the number of matches within 3" = 0.0008degrees
   index = where(dist lt 0.0008, count)
   sep=min(dist,ind)
    if (sep LT 0.0003)  then begin
      mmatch[q]=ind

      printf, outlun2, count, magauto(q) - object[ind].imaga
      countarr[cind] = count
      deltaarr[cind] = magauto(q) - object[ind].imaga
      cind = cind + 1
   endif
endfor

print, countarr, deltaarr
matched=where(mmatch GT 0)
nonmatched = where(mmatch lt 0)

print, n_elements(matched),"matched"
print, n_elements(nonmatched),"nonmatched"
z = 0
for num = 0, n_elements(matched) - 1 do begin

   iarr(z) = object(mmatch[matched[num]]).imaga
   magisoarr(z)= magiso(matched[num])
   magisocorarr(z)= magisocor(matched[num])
   magaperarr(z)= magaper(matched[num])
   magautoarr(z)= magauto(matched[num])
   magbestarr(z)= magbest(matched[num])
   fwhmarr(z) = fwhm(matched[num])
   xarr(z) = x(matched[num])
   yarr(z) = y(matched[num])
   raarr(z) = ra(matched[num])
   decarr(z) = dec(matched[num])

   z = z + 1
endfor

plot, iarr, magautoarr, psym = 2, xrange=[10,30], yrange=[10,30], xtitle='palomari', ytitle='magauto'
plot, iarr, magautoarr - iarr, psym = 2, xrange=[10,30], yrange=[-5,5], xtitle='palomari', ytitle='magauto - palomar'

plot, iarr, magisoarr, psym = 2, xrange=[10,30], yrange=[10,30], xtitle='palomari', ytitle='magiso'
plot, iarr, magisoarr - iarr, psym = 2, xrange=[10,30], yrange=[-5,5], xtitle='palomari', ytitle='magiso - palomar'

plot, iarr, magisocorarr, psym = 2, xrange=[10,30], yrange=[10,30], xtitle='palomari', ytitle='magisocor'
plot, iarr, magisocorarr - iarr, psym = 2, xrange=[10,30], yrange=[-5,5], xtitle='palomari', ytitle='magisocor - palomar'

plot, iarr, magaperarr, psym = 2, xrange=[10,30], yrange=[10,30], xtitle='palomari', ytitle='magaper'
plot, iarr, magaperarr - iarr, psym = 2, xrange=[10,30], yrange=[-5,5], xtitle='palomari', ytitle='magaper - palomar'

plot, iarr, magbestarr, psym = 2, xrange=[10,30], yrange=[10,30], xtitle='palomari', ytitle='magbest'
plot, iarr, magbestarr - iarr, psym = 2, xrange=[10,30], yrange=[-5,5], xtitle='palomari', ytitle='magbest - palomar'

plot, magauto, magbest, psym = 2, xrange=[10,30], yrange=[10,30], xtitle='magauto', ytitle='magbest'

plot, deltaarr, countarr, psym = 2, xrange=[-5,5], thick = 3
plothist, countarr

openw, outlun, '/Users/jkrick/hst/xy.out', /get_lun

goodarr = fltarr(n_elements(iarr))
badarr = fltarr(n_elements(iarr))
k = 0
l = 0
for j = 0, n_elements(iarr) - 1 do begin
   if abs(magbestarr(j) - iarr(j)) lt 1.0 then begin
      goodarr(k) = j
      k = k + 1
   endif else begin
      badarr(l) = j
      printf, outlun, raarr(j), decarr(j)
      l = l + 1
   endelse

endfor
goodarr = goodarr[0:k-1]
badarr =badarr[0:l-1]



ps_close, /noprint,/noid
close, outlun
free_lun, outlun
close, outlun2
free_lun, outlun2

end
