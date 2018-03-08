pro iracphot_check
!p.multi = [0, 2, 2]

restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'
nflux = get_nflux()

;want to look at the 'isolated galaxies'
;can use nearest neighbor distance as a discriminator
a = where(objectnew.irac1flux gt 0 and objectnew.irac2flux gt 0 and objectnew.irac3flux gt 0 and objectnew.irac4flux gt 0 and objectnew.nearestdist *3600. gt 4.); and nflux gt 5)
;plothist, objectnew[a].nearestdist *3600.
print, 'a', n_elements(a)

plot, objectnew[a].irac1flux, objectnew[a].irac1fluxold, psym = 2, xrange=[0,50], yrange=[0,50], xtitle = 'photometry with ACS priors', ytitle = 'IRAC only photometry'
oplot, findgen(100)*1E3, findgen(100)*1E3 

plot, objectnew[a].irac2flux, objectnew[a].irac2fluxold, psym = 2, xrange=[0,50], yrange=[0,50], xtitle = 'photometry with ACS priors', ytitle = 'IRAC only photometry'
oplot, findgen(100)*1E3, findgen(100)*1E3

plot, objectnew[a].irac3flux, objectnew[a].irac3fluxold, psym = 2, xrange=[0,50], yrange=[0,50], xtitle = 'photometry with ACS priors', ytitle = 'IRAC only photometry'
oplot, findgen(100)*1E3, findgen(100)*1E3

plot, objectnew[a].irac4flux, objectnew[a].irac4fluxold, psym = 2, xrange=[0,50], yrange=[0,50], xtitle = 'photometry with ACS priors', ytitle = 'IRAC only photometry'
oplot, findgen(100)*1E3, findgen(100)*1E3

ps_open, filename='/Users/jkrick/nutella/nep/datapaper/photerr_laidler.ps',/portrait,/square,/color

!p.multi = [0, 0, 1]
vsym, /polygon, /fill

a = where(nflux gt 5); objectnew.nearestdist *3600. gt 0. and nflux gt 5)
print, 'n_a', n_elements(a)
;plot, objectnew[a].irac1magold, (objectnew[a].irac1magold - objectnew[a].irac1mag) , psym = 8, xtitle='IRAC only magnitudes (AB)', ytitle='IRAC only - ACS prior derived mags',  charthick=3, xthick=3, ythick=3, symsize = 0.3

;plot, objectnew[a].irac1fluxold, (objectnew[a].irac1flux - objectnew[a].irac1fluxold)/objectnew[a].irac1fluxold , psym = 8, xtitle='IRAC only flux', ytitle='(ACS prior derived - IRAC only) / irac only',  charthick=3, xthick=3, ythick=3, symsize = 0.3, xrange= [1, 2E4], /xlog, yrange = [-5,5]

;oplot, objectnew[a].irac1fluxold, fltarr(n_elements(a))
j = where(objectnew.acsmag lt 0 or objectnew.acsmag gt 90)
i = where(objectnew[j].irac1mag gt 0 and objectnew[j].irac1mag lt 90)


;want to bin these by irac mag.
;error bars by the stddev within the bin.
x = objectnew[a].irac1fluxold
y =  (objectnew[a].irac1flux - objectnew[a].irac1fluxold)/objectnew[a].irac1fluxold 

binsize = 10
xarr = fltarr(5000)
yarr = fltarr(5000)
yerrarr = fltarr(5000)
count = 0
for n = 1., 100, binsize do begin
   print, 'n ', n
   bin1 = where(x gt n and x le n+binsize)
;   print, 'x(bin1)', x(bin1)
;   print, 'y(bin1)', y(bin1)
   print, '(n + n+binsize) / 2', (n + n+binsize) / 2
   print, n_elements(bin1)
   xarr(count) = (n + n+binsize) / 2
   mean_biweight = biweight_mean(y(bin1), sigma_biweight)
   yarr(count) = mean_biweight; mean(y(bin1))
   yerrarr(count) = sigma_biweight; sigma_robust; stddev(y(bin1))
   count = count + 1
endfor

binsize = 100
for n = 100., 1000, binsize do begin
   print, 'n ', n
   bin1 = where(x gt n and x le n+binsize)
;   print, 'x(bin1)', x(bin1)
;   print, 'y(bin1)', y(bin1)
   print, '(n + n+binsize) / 2', (n + n+binsize) / 2
   print, n_elements(bin1)
   mean_biweight = biweight_mean(y(bin1), sigma_biweight)
   xarr(count) = (n + n+binsize) / 2
   yarr(count) = mean_biweight              ;mean(y(bin1))
   yerrarr(count) = sigma_biweight ; sigma_robust; robust_sigma(y(bin1))
   count = count + 1
endfor

xarr = xarr[0:count - 1]
yarr = yarr[0:count - 1]
yerrarr = yerrarr[0:count -1]
xerrarr = fltarr(n_elements(yerrarr))
print, 'yarr',n_elements(yerrarr)
xerrarr[0:9] = 5
xerrarr[10:19] = 50
print, 'err', xerrarr 
plot, xarr, yarr, yrange = [-1,1], xrange = [1, 1E3], /xlog,ystyle = 1, psym =2, $
      xtitle = 'F(IRAC) in microJy', ytitle = '(F(ACS) - F(IRAC)) / F(IRAC)', charthick = 3, xthick = 3, ythick = 3, thick = 3
;oploterr, xarr, yarr,  yerrarr
oploterror, xarr, yarr, xerrarr, yerrarr, psym = 2

;oplot, findgen(1E5), fltarr(1E5)
;irac extended
ei1 = where(objectnew[j[i]].extended gt 0)
pi1 = where(objectnew[j[i]].extended lt 1)
;oplot, objectnew[j[i[ei1]]].irac1magold, objectnew[j[i[ei1]]].irac1magold - objectnew[j[i[ei1]]].irac1mag, psym = 2
;oplot, objectnew[j[i[pi1]]].irac1magold, objectnew[j[i[pi1]]].irac1magold - objectnew[j[i[pi1]]].irac1mag, psym = 4


a = where(objectnew[a].irac1magold - objectnew[a].irac1mag gt 0.5 and objectnew.irac1magold gt 24 and objectnew.irac1magold lt 90)

;openw, outlunred, '/Users/jkrick/nep/datapaper/iracphot.reg', /get_lun
;printf, outlunred, 'fk5'
;for rc=0, n_elements(a) -1 do  printf, outlunred, 'circle( ', objectnew[a[rc]].ra, objectnew[a[rc]].dec, ' 3")'
;close, outlunred
;free_lun, outlunred

ps_close, /noprint,/noid

end
