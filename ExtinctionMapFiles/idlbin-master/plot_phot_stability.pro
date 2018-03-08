pro plot_phot_stability, chname
dirname =  '/Users/jkrick/irac_warm/calstars/pmap_star_' + chname + '/'
if chname eq 'ch2' then starname = 'bd67_1044' else starname = 'KF09T1'

savename = strcompress(dirname + starname + '_phot_'+chname+'.sav',/remove_all)
print, 'savename', savename
restore, savename

cd,dirname
readcol, 'aorlist.txt', aorname, format = 'A', /silent
if chname eq 'ch2' then normflux = 0.399 else normflux = 0.146
if chname eq 'ch2' then time0 = 5.57E4 else time0 = 5.58E4

bmjdg1arr = fltarr(n_elements(aorname))
fluxg1arr = bmjdg1arr
fluxerrg1arr = bmjdg1arr
xarrg1arr = bmjdg1arr
yarrg1arr = bmjdg1arr
bmjdg2arr = bmjdg1arr
fluxg2arr = bmjdg1arr
fluxerrg2arr = bmjdg1arr
xarrg2arr = bmjdg1arr
yarrg2arr = bmjdg1arr
bmjdg3arr = bmjdg1arr
fluxg3arr = bmjdg1arr
fluxerrg3arr = bmjdg1arr
xarrg3arr = bmjdg1arr
yarrg3arr = bmjdg1arr
bmjdg4arr = bmjdg1arr
fluxg4arr = bmjdg1arr
fluxerrg4arr = bmjdg1arr
xarrg4arr = bmjdg1arr
yarrg4arr = bmjdg1arr

totalg1count = 0
totalg2count = 0
totalg3count = 0
totalg4count = 0
c1 = 0
c2 = 0
c3 = 0
c4 = 0

nimages = 64  ; 10

if chname eq 'ch2' then fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/occu_ch2_0p1s_x4_500x500_0043_121120.fits', pmapdata, pmapheader
if chname eq 'ch1' then fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/occu_ch1_500x500_0043_120409.fits', pmapdata, pmapheader
c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', title = planetname, aspect_ratio = 1, xrange = [0,500], yrange = [0,500])

for a = 0, n_elements(aorname) -1 do begin
   print, 'working on aorname ', aorname(a) , a, n_elements(aorname)

   ;sort out the straight up bad data points.
   bmjd= planethash[aorname(a),'bmjdarr']
   flux = planethash[aorname(a),'fluxarr']
   fluxerr = planethash[aorname(a),'fluxerrarr']
   xarr = planethash[aorname(a),'xarr']
   yarr = planethash[aorname(a),'yarr']

   good = where (finite(flux) gt 0 and bmjd gt 0,complement = nbad)
  ; print, 'bad', n_elements(nbad), n_elements(flux)
 ;  if n_elements(nbad) gt 1 then print, 'testing bad', flux(nbad(0)), flux(nbad(1))
   bmjd = bmjd(good)
   flux = flux(good)
   fluxerr = fluxerr(good)
   xarr = xarr(good)
   yarr = yarr(good)


   ;select only a narrow range in position
   if chname eq 'ch1' then begin
      ;looking for more counts
      g1 = where(xarr ge 15.1 and xarr le 15.14 and yarr ge 15.1 and yarr le 15.14, g1count)   ;black
      g2 = where(xarr ge 15.14 and xarr le 15.18 and yarr ge 15.1 and yarr le 15.14, g2count)  ;red
      g3 = where(xarr ge 15.10 and xarr le 15.14 and yarr ge 15.06 and yarr le 15.1, g3count)  ;blue
      g4 = where(xarr ge 15.14 and xarr le 15.18 and yarr ge 15.06 and yarr le 15.1, g4count) ;cyan
   endif else begin ; ch2
      ;smaller
      g1 = where(xarr ge 15.1 and xarr le 15.12 and yarr ge 15.1 and yarr le 15.12, g1count) ;black
      g2 = where(xarr ge 15.13 and xarr le 15.15 and yarr ge 15.1 and yarr le 15.12, g2count) ;red
      g3 = where(xarr ge 15.10 and xarr le 15.12 and yarr ge 15.07 and yarr le 15.09, g3count) ;blue
      g4 = where(xarr ge 15.13 and xarr le 15.15 and yarr ge 15.07 and yarr le 15.09, g4count) ;cyan
 ;    g1 = where(xarr ge 15.1 and xarr le 15.14 and yarr ge 15.1 and yarr le 15.14, g1count) ;black
 ;     g2 = where(xarr ge 15.14 and xarr le 15.18 and yarr ge 15.1 and yarr le 15.14, g2count) ;red
 ;     g3 = where(xarr ge 15.10 and xarr le 15.14 and yarr ge 15.06 and yarr le 15.1, g3count) ;blue
 ;     g4 = where(xarr ge 15.14 and xarr le 15.18 and yarr ge 15.06 and yarr le 15.1, g4count) ;cyan

   endelse



;   print, 'narrow position', g1count
   if g1count gt nimages then begin
      totalg1count = totalg1count + g1count
      bmjdg1arr[c1] = median(bmjd(g1)) 
      fluxg1arr[c1] = mean(flux(g1)/ normflux)
      fluxerrg1arr[c1] = stddev(flux(g1)/ normflux)
      xarrg1arr[c1] = mean(xarr(g1))
      yarrg1arr[c1] = mean(yarr(g1))
      xcen = 500.* (xarr(g1)- 14.5)
      ycen = 500.* (yarr(g1) - 14.5)
      pc1 = plot(xcen, ycen, '1s', sym_size = 0.1,   sym_filled = 1,  overplot = c, color = 'black')
      c1++

   endif
   if g2count gt nimages then begin
      totalg2count = totalg2count + g2count
      bmjdg2arr[c2] = median(bmjd(g2)) 
      fluxg2arr[c2] = mean(flux(g2)/ normflux)
      fluxerrg2arr[c2] = stddev(flux(g2)/ normflux)
      xarrg2arr[c2] = mean(xarr(g2))
      yarrg2arr[c2] = mean(yarr(g2))
      xcen = 500.* (xarr(g2)- 14.5)
      ycen = 500.* (yarr(g2) - 14.5)
      pc2 = plot(xcen, ycen, '1s', sym_size = 0.1,   sym_filled = 1,  overplot = c, color = 'red')
      c2++
   endif
  if g3count gt nimages then begin
     totalg3count = totalg3count + g3count
      bmjdg3arr[c3] = median(bmjd(g3)) 
      fluxg3arr[c3] = mean(flux(g3)/ normflux)
      fluxerrg3arr[c3] = stddev(flux(g3)/ normflux)
      xarrg3arr[c3] = mean(xarr(g3))
      yarrg3arr[c3] = mean(yarr(g3))
      xcen = 500.* (xarr(g3)- 14.5)
      ycen = 500.* (yarr(g3) - 14.5)
      pc3 = plot(xcen, ycen, '1s', sym_size = 0.1,   sym_filled = 1,  overplot = c, color = 'blue')
      c3++
   endif
  if g4count gt nimages then begin
     totalg4count = totalg4count + g4count
      bmjdg4arr[c4] = median(bmjd(g4)) 
      fluxg4arr[c4] = mean(flux(g4)/ normflux)
      fluxerrg4arr[c4] = stddev(flux(g4)/ normflux)
      xarrg4arr[c4] = mean(xarr(g4))
      yarrg4arr[c4] = mean(yarr(g4))
      xcen = 500.* (xarr(g4)- 14.5)
      ycen = 500.* (yarr(g4) - 14.5)
      pc4 = plot(xcen, ycen, '1s', sym_size = 0.1,   sym_filled = 1,  overplot = c, color = 'cyan')
      c4++
   endif

endfor

bmjdg1arr = bmjdg1arr[0:c1 - 1]
fluxg1arr = fluxg1arr[0:c1-1] 
fluxerrg1arr= fluxerrg1arr[0:c1 - 1]
xarrg1arr = xarrg1arr[0:c1-1] 
yarrg1arr = yarrg1arr[0:c1-1] 

bmjdg2arr = bmjdg2arr[0:c2 - 1]
fluxg2arr = fluxg2arr[0:c2-1] 
fluxerrg2arr= fluxerrg2arr[0:c2 - 1]

bmjdg3arr = bmjdg3arr[0:c3 - 1]
fluxg3arr = fluxg3arr[0:c3-1] 
fluxerrg3arr= fluxerrg3arr[0:c3 - 1]

bmjdg4arr = bmjdg4arr[0:c4 - 1]
fluxg4arr = fluxg4arr[0:c4-1] 
fluxerrg4arr= fluxerrg4arr[0:c4 - 1]

print, 'counts', totalg1count, totalg2count, totalg3count, totalg4count

;------------------plot it
print, 'testing ', fluxg1arr
p1 = plot(bmjdg1arr - time0,fluxg1arr,  title = starname + chname, xtitle = 'Time(days)', ytitle = 'Normalized Mean Flux', '1s', sym_size = 0.5,   sym_filled = 1, xthick = 2, ythick = 2) ;, errorbar_capsize = 0.005, , xrange = [0,600], yrange = [0.994, 1.006]
;t1 = text(450, 0.998, string(totalg1count),/data)

p2 = plot(bmjdg2arr - time0,fluxg2arr,   '1s', sym_size = 0.5,   sym_filled = 1, overplot = p1, color = 'red' )
;t2 = text(450, 0.996, string(totalg2count),/data, color = 'red')

p3 = plot(bmjdg3arr - time0,fluxg3arr,   '1s', sym_size = 0.5,   sym_filled = 1, overplot = p1, color = 'blue' )
;t3 = text(450, 0.994, string(totalg3count),/data, color = 'blue')

p4 = plot(bmjdg4arr - time0,fluxg4arr,   '1s', sym_size = 0.5,   sym_filled = 1, overplot = p1, color = 'cyan' )
;t4 = text(450, 0.992, string(totalg4count),/data, color = 'cyan')
;----------------------------------------------
;fit with a liner function.

start = [1E-5,1.0]
startflat = [1.0]

result= MPFITFUN('linear',bmjdg1arr - time0,fluxg1arr, fluxerrg1arr, start, perror = perror, bestnorm = bestnorm)    
fp = plot(bmjdg1arr - time0, result(0)*(bmjdg1arr - time0) + result(1), $
   color = 'black', overplot = p1,thick = 3)
DOF     = N_ELEMENTS(bmjdg1arr) - N_ELEMENTS(result) ; deg of freedom
PCERROR = PERROR * SQRT(BESTNORM / DOF)     ; scaled uncertainties
print, 'black result', result, 'pcerror', pcerror, ' perror', perror, 'bestnorm', bestnorm, 'dof', DOF, sqrt(bestnorm/dof)
t1 = text(150, 0.997, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm, 5), /data, color = 'black', font_style = 'bold');pcerror(0)

result2 = MPFITFUN('flatline',bmjdg1arr - time0,fluxg1arr, fluxerrg1arr, startflat, perror = perror, bestnorm = bestnorm)    
;fp = plot(bmjdg1arr - time0, 0*(bmjdg1arr - time0) + result2(0), $
;   color = 'black', overplot = p1,thick = 3, linestyle = 2)
;tt1 = text(400, 0.997, '0    ' + sigfig(bestnorm, 5), /data, color = 'black', font_style = 'bold')
;--------------
result= MPFITFUN('linear',bmjdg2arr - time0,fluxg2arr, fluxerrg2arr, start, perror = perror, bestnorm = bestnorm,/quiet)    
fp = plot(bmjdg2arr - time0, result(0)*(bmjdg2arr - time0) + result(1), $
   thick = 3, color = 'red', overplot = p1)
DOF     = N_ELEMENTS(bmjdg2arr) - N_ELEMENTS(result) ; deg of freedom
PCERROR = PERROR * SQRT(BESTNORM / DOF)     ; scaled uncertainties
print, '-------------------------'
print, 'red result', result, 'pcerror', pcerror, ' perror', perror
t1 = text(150, 0.996, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm, 5), /data, color = 'red', font_style = 'bold')

result2 = MPFITFUN('flatline',bmjdg2arr - time0,fluxg2arr, fluxerrg2arr, startflat, perror = perror, bestnorm = bestnorm)    
;fp = plot(bmjdg2arr - time0, 0*(bmjdg2arr - time0) + result2(0), $
;   color = 'red', overplot = p1,thick = 3, linestyle = 2)
;tt1 = text(400, 0.9965, '0    ' + sigfig(bestnorm, 5), /data, color = 'red', font_style = 'bold')
;--------------


result= MPFITFUN('linear',bmjdg3arr - time0,fluxg3arr, fluxerrg3arr, start, perror = perror, bestnorm = bestnorm,/quiet)    
fp = plot(bmjdg3arr - time0, result(0)*(bmjdg3arr - time0) + result(1), $
   thick = 3, color = 'blue', overplot = p1)
DOF     = N_ELEMENTS(bmjdg3arr) - N_ELEMENTS(result) ; deg of freedom
PCERROR = PERROR * SQRT(BESTNORM / DOF)     ; scaled uncertainties
print, '-------------------------'
print, 'blue result', result, 'pcerror', pcerror, ' perror', perror
t1 = text(150, 0.995, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm, 5), /data, color = 'blue', font_style = 'bold')

result2 = MPFITFUN('flatline',bmjdg3arr - time0,fluxg3arr, fluxerrg3arr, startflat, perror = perror, bestnorm = bestnorm)    
;fp = plot(bmjdg3arr - time0, 0*(bmjdg3arr - time0) + result2(0), $
;   color = 'blue', overplot = p1,thick = 3, linestyle = 2)
;tt1 = text(400, 0.996, '0    ' + sigfig(bestnorm, 5), /data, color = 'blue', font_style = 'bold')
;--------------

result= MPFITFUN('linear',bmjdg4arr - time0,fluxg4arr, fluxerrg4arr, start, perror = perror, bestnorm = bestnorm,/quiet)   
fp = plot(bmjdg4arr - time0, result(0)*(bmjdg4arr - time0) + result(1), $
   thick =3, color = 'cyan', overplot = p1)
DOF     = N_ELEMENTS(bmjdg4arr) - N_ELEMENTS(result) ; deg of freedom
PCERROR = PERROR * SQRT(BESTNORM / DOF)     ; scaled uncertainties
print, '-------------------------'
print, 'cyan result', result, 'pcerror', pcerror, ' perror', perror
t1 = text(150, 0.994, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm, 5), /data, color = 'cyan', font_style = 'bold')

result2 = MPFITFUN('flatline',bmjdg4arr - time0,fluxg4arr, fluxerrg4arr, startflat, perror = perror, bestnorm = bestnorm)    
;fp = plot(bmjdg4arr - time0, 0*(bmjdg4arr - time0) + result2(0), $
;   color = 'cyan', overplot = p1,thick = 3, linestyle = 2)
;tt1 = text(400, 0.9955, '0    ' + sigfig(bestnorm, 5), /data, color = 'cyan', font_style = 'bold')
;--------------

end
