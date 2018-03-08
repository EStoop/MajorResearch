pro plot_wasp14_ch2
  colorarr = ['blue', 'red', 'aquamarine',  'dark_blue', 'dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']

  
  restore,'/Users/jkrick/irac_warm/pcrs_planets/wasp14/wasp14_phot_ch2.sav'
  aorname = ['r45426688',  'r45428224',  'r45428480',  'r45428736',  'r45428992'] ;ch2


;-------------------------------------------------------
;taken from nsted planet properties calculator
  utmjd_center = double(56086.74059); 2456087.24059

  duration = 226.4             ;transit duration in minutes
  duration = duration /60./24.  ; in days
  period = 2.519961 ;days
;set phase = 0.0 in the middle of the transit
;  for a = 0, n_elements(aorname) - 1 do  begin
;      AORwasp14[a].bmjdarr = ((AORwasp14[a].bmjdarr - utmjd_center) / period) 
;   endfor

;print, 'phase', AORwasp14[0].bmjdarr[0:100]


;  fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader
;  c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500])
;  for a = 0, n_elements(aorname) - 1 do begin
;     xcen500 = 500.* (AORwasp14[a].xcen - 14.5)
;     ycen500 = 500.* (AORwasp14[a].ycen - 14.5)
;     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
; endfor

;-----
 
;    for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'X pix',  title = 'wasp14 ch2', xrange = [0,60], yrange = [14.5, 15.5])
;     endif else begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------

;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Y pix',  title = 'wasp14 ch2', xrange = [0,60], yrange = [14.5, 15.5])
;     endif else begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Background', xrange = [0,10])
;     endif else begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].flux, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Flux',  title = 'wasp14 ch2', xrange = [0,60], yrange = [0.057, 0.063])
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].corrflux-0.002,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)
;     endif else begin
;        am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot)
; am = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].corrflux-0.002,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)
;     endelse

;  endfor

;------
;a = where(finite(AORwasp14[1].corrflux) gt 0 and finite(AORwasp14[1].flux) gt 0)

; print, (abs(AORwasp14[1].corrflux[a[0:100]] - AORwasp14[1].flux[a[0:100]])/ AORwasp14[1].flux[a[0:100]])*100
;plothist,(abs(AORwasp14[1].corrflux[a[0:100]] - AORwasp14[1].flux[a[0:100]])/ AORwasp14[1].flux[a[0:100]])*100, xhist, yhist, bin = 0.1,/fill


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;try this diffferently, try putting the AORs togetether as one array; then bin

timearr = [AORwasp14[0].timearr, AORwasp14[1].timearr, AORwasp14[2].timearr, AORwasp14[3].timearr, AORwasp14[4].timearr]
fluxarr = [AORwasp14[0].flux, AORwasp14[1].flux, AORwasp14[2].flux, AORwasp14[3].flux, AORwasp14[4].flux]
fluxerrarr = [AORwasp14[0].fluxerr, AORwasp14[1].fluxerr, AORwasp14[2].fluxerr, AORwasp14[3].fluxerr, AORwasp14[4].fluxerr]
corrfluxarr = [AORwasp14[0].corrflux, AORwasp14[1].corrflux, AORwasp14[2].corrflux, AORwasp14[3].corrflux, AORwasp14[4].corrflux]
corrfluxerrarr = [AORwasp14[0].corrfluxerr, AORwasp14[1].corrfluxerr, AORwasp14[2].corrfluxerr, AORwasp14[3].corrfluxerr, AORwasp14[4].corrfluxerr]
xarr = [AORwasp14[0].xcen, AORwasp14[1].xcen, AORwasp14[2].xcen, AORwasp14[3].xcen, AORwasp14[4].xcen]
yarr = [AORwasp14[0].ycen, AORwasp14[1].ycen, AORwasp14[2].ycen, AORwasp14[3].ycen, AORwasp14[4].ycen]
bkgd = [AORwasp14[0].bkgd, AORwasp14[1].bkgd, AORwasp14[2].bkgd, AORwasp14[3].bkgd, AORwasp14[4].bkgd]
bmjd = [AORwasp14[0].bmjdarr, AORwasp14[1].bmjdarr, AORwasp14[2].bmjdarr, AORwasp14[3].bmjdarr, AORwasp14[4].bmjdarr]

print, 'n', n_elements(timearr), n_elements(AORwasp14[0].timearr), n_elements(AORwasp14[1].timearr)
;now for the binning

      ;remove outliers
good = where(xarr lt mean(xarr) + 5.*stddev(xarr) and xarr gt mean(xarr) -5.*stddev(xarr) and xarr lt mean(xarr) +5.0*stddev(yarr) and yarr gt mean(yarr) - 5.0*stddev(yarr),ngood_pmap, complement=bad) 
     print, 'bad ',n_elements(bad)
     xarr = xarr[good]
     yarr = yarr[good]
     timearr = timearr[good]
     flux = fluxarr[good]
     fluxerr = fluxerrarr[good]
     corrflux = corrfluxarr[good]
     corrfluxerr = corrfluxerrarr[good]
     bmjdarr = bmjd[good]
     bkgdarr = bkgd[good]
     
     
; binning
     bin_level = 63*1L;*60
     nframes = bin_level
     nfits = long(n_elements(corrflux)) / nframes
     print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
     bin_flux = dblarr(nfits)
     bin_fluxerr = dblarr(nfits)
     bin_corrflux= dblarr(nfits)
     bin_corrfluxerr= dblarr(nfits)
     bin_timearr = dblarr(nfits)
     bin_bmjdarr = dblarr(nfits)
     bin_bkgd = dblarr(nfits)

     for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', n_elements(AORwasp14[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
        idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
        id = flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
        ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
       ; if si lt 2 then print, 'ib', ib
       ; print, 'mean', mean(ib,/nan)
        bin_flux[si] = mean(id, /nan)
        bin_fluxerr[si] = sqrt(total(id^2))/ (n_elements(id))
        bin_corrflux[si] = mean(idata,/nan)
        bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        bin_timearr[si]= mean(timearr[si*nframes:si*nframes + (nframes - 1)])
        bin_bmjdarr[si]= mean(bmjdarr[si*nframes:si*nframes + (nframes - 1)])
        bin_bkgd[si] = mean(ib,/nan)

     endfor                     ;for each fits image
     
 ;try getting rid of flux outliers.
  ;do some running mean with clipping
     start = 0
     print, 'nflux', n_elements(bin_flux)
     for ni = 50, n_elements(bin_flux) -1, 50 do begin
        meanclip,flux[start:ni], m, s, subs = subs ;,/verbose
                                ;print, 'good', subs+start
                                ;now keep a list of the good ones
        if ni eq 50 then good_ni = subs+start else good_ni = [good_ni, subs+start]
        start = ni
     endfor
                                ;see if it worked
     ;xarr = xarr[good_ni]
     ;yarr = yarr[good_ni]
     bin_timearr = bin_timearr[good_ni]
     bin_bmjdarr = bin_bmjdarr[good_ni]
     bin_corrfluxerr = bin_corrfluxerr[good_ni]
     bin_flux = bin_flux[good_ni]
     bin_fluxerr = bin_fluxerr[good_ni]
     bin_corrflux = bin_corrflux[good_ni]
     bin_bkgd = bin_bkgd[good_ni]

     print, 'nflux', n_elements(bin_flux), n_elements(bin_timearr)
     
save, bin_timearr, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_timearr.sav'
save, bin_flux, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_flux.sav'
save, bin_corrflux, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_corrflux.sav'

;now try plotting
;     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux/median(bin_flux), '6r1s', sym_size = 0.1,   sym_filled = 1,  yrange = [0.975, 1.01], xtitle = 'Time (hrs)', ytitle = 'Normalized Binned Flux',  title = 'wasp14 ch2', xrange = [0,60])
;     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux/median(bin_corrflux)-0.01, '61s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)

;     pl = errorplot(bin_bmjdarr, bin_flux/median(bin_flux) ,bin_fluxerr, '6r1s', sym_size = 0.1,   sym_filled = 1,   xtitle = 'Phase', ytitle = 'Normalized Flux', errorbar_capsize = 0., xrange = [-1,1], yrange = [0.99, 1.01])
;     pl = errorplot(bin_bmjdarr, (bin_corrflux/median(bin_corrflux))-0.003, bin_corrfluxerr,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot, errorbar_capsize = 0.)

;     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_bkgd/median(bin_bkgd), '61s', sym_size = 0.1,   sym_filled = 1,  xtitle = 'Time (hrs)', ytitle = 'Normalized binned Bkgd',  title = 'wasp14 ch2', xrange = [0,60])


;------------------------------------------------------------
;periodogram
;print, 'starting periodogram'
;X = (bin_timearr - bin_timearr(0))/60./60. ;in hours
X = (bin_timearr - bin_timearr(0))/60./60. ;in days
Y = bin_corrflux
;t = plot(x, y)
;ok, try cutting this down in phase to between eclipses
;a = where(x gt 7 and x lt 17)  ;in hours
a1 = where(x gt 7.5 and x lt 31.) ; in hours
a2 = where(x gt 35.5 and x lt 44.) ; in days 0.9

x1 =x(a1)
x2 = x(a2)
y1 = y(a1)
y2=y(a2)
ytot = [y1, y2]

;
;try taking out the phase variation
;which also takes out the mean
start =  [-0.001,0.173]
noise1 = fltarr(n_elements(x1)) + 1. ;set all noises equal
print, 'n fit', n_elements(x1), n_elements(y1), n_elements(noise1)
fitres1= MPFITFUN('linear',x1,y1, noise1, start)    ;ICL
print, 'fitres1', fitres1
noise2 = fltarr(n_elements(x2)) + 1. ;set all noises equal
fitres2= MPFITFUN('linear',x2,y2, noise2, start)    ;ICL
ap = plot(x1, y1)
ap = plot(x1, fitres1(0)*x1 + fitres1(1),color = 'green',/overplot)
a = plot(x2, y2, /overplot)
ap = plot(x2, fitres2(0)*x2 + fitres2(1),color = 'green',/overplot)
;
flaty1 = y1 - (fitres1(0)*x1 + fitres1(1))
flaty2 = y2 - (fitres2(0)*x2+ fitres2(1))
x = [x1, x2]
flaty = [flaty1, flaty2]
ap = plot(x, flaty+0.06, color = 'blue', /overplot)

;now put them together after flattening
;print, y[0:100]
; Test the hypothesis that X and Y represent a significant periodic
; signal against the hypothesis that they represent random noise:
result = LNP_TEST(X, flaty,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
;PRINT, 'lnp result', result

b = plot(wk1, wk2, xtitle = 'Frequency(1/days)', ytitle = 'Power', thick = 2, name = 'PMAP Corrected', xrange = [0,100]);

;------------------------------------------------------------
;now try overplotting uncorrected periodogram
Y2 = bin_flux
;ok, try cutting this down in phase to between eclipses
y2 = [Y2(a1),Y2(a2)]
badx = where(finite(x) lt 1, badxcount)
bady = where(finite(y2) lt 1, badycount)
print,'badycount', badycount
if badycount gt 0 then y2(bady) = y2(bady-1) 
;print, y2[0:100]

;set the mean to zero
;doesn't seem to make a difference.
;y2 = y2 - mean(y2)
;print, y2[0:100]

;try taking out the phase variation
;which also takes out the mean
start = [0.001,0.001]
noise = fltarr(n_elements(x)) + 1. ;set all noises equal
fitres= MPFITFUN('linear',x,y2, noise, start)    ;ICL
;ap = plot(x, y2)
;ap = plot(x, fitres(0)*x + fitres(1),color = 'green',/overplot)

flaty2 = y2- (fitres(0)*x + fitres(1))
;ap = plot(x, flaty2, color = 'blue', /overplot)

; Test the hypothesis that X and Y represent a significant periodic
; signal against the hypothesis that they represent random noise:
result = LNP_TEST(X, flatY2,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
PRINT, result

b2 = plot(wk1, wk2, /overplot, color = 'red',thick = 2, name = 'Raw')
;---------------------------------------------------------------------
;XXX just need to do this 1000 times now.

;now try monte carloing the noise periodogram
;randomize the order of the photometry while keeping the time
;information the same.
nwk2 = n_elements(wk2) ;which is bin_flux cut down
ny2= n_elements(y2)
;make a randomly ordered array with nel elements
nmc = 100
randresult = fltarr(nmc, nwk2)
for rcount = 0, nmc - 1 do begin
;   print, 'rcount', rcount
   rand = randomu(seed, ny2)
;now use that as the pointer/orderer for the photometry
   randflux = y2[sort(rand)]
   rr = LNP_TEST(X, randflux, /double,  WK1 = wk1, WK2 = wk2, JMAX = jmax)
;   print, 'n wk2', n_elements(wk2), n_elements(x), n_elements(y2)
   randresult[rcount,*] = wk2
  ; randp = plot(wk1, wk2,/overplot, color = 'magenta')

;   randflux = ytot[sort(rand)]
;   rr  = LNP_TEST(X, randflux, /double,  WK1 = wk1, WK2 = wk2, JMAX = jmax)
;   randresult[rcount,*] = wk2
   ;randp = plot(wk1, wk2,/overplot, color = 'grey')
endfor
final_max = fltarr(nwk2)
final_mean = fltarr(nwk2)
final_min = fltarr(nwk2)

for rfreq = 0, nwk2 -1 do begin
  ; print, 'rfreq', rfreq
   ;help, a
   final_max(rfreq) = max(randresult[*,rfreq])
   final_mean(rfreq) = mean(randresult[*,rfreq])
   final_min(rfreq) = min(randresult[*,rfreq])
endfor
;cut down to fit on the plot
good = where(wk1 le 50)
;d = plot(wk1(good), final_max(good), /overplot, color = 'black')
;d = plot(wk1(good), final_mean(good), /overplot, color = 'black')
;d = plot(wk1, final_min, /overplot, color = 'cyan')

;poly = polygon([wk1(good), reverse(wk1(good))], [final_max(good),reverse(final_mean(good))], /data, /fill_background, fill_color ='grey', name = 'noise')
;l = legend(target = [b2, b, poly], position = [35, 140],/data)

;------------------------------------------------------------
; and the periodogram for the self calibrated data

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_time.sav'
restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_flux.sav'
x = bin_time /24.
y = bin_flux

a1 = where(x gt 7.1/24 and x lt 30.5/24.) ; in days
a2 = where(x gt 35.5/24. and x lt 44./24.)
x1 =x(a1)
x2 = x(a2)
y1 = y(a1)
y2=y(a2)

;badx = where(finite(x) lt 1, badxcount)
;bady = where(finite(y1) lt 1, badycount)
;print,'badycount', badycount
;print, x(bady), y(bady)
;if badycount gt 0 then y1(bady) = y1(bady-1) 
;print, x(bady), y(bady)
;print, y[0:100]
;ytot = [y1, y2,y3,y4,y5]

;try taking out the phase variation
;which also takes out the mean
start =  [-0.001,0.173]
noise1 = fltarr(n_elements(x1)) + 1. ;set all noises equal
fitres1= MPFITFUN('linear',x1,y1, noise1, start)    ;ICL
print, 'fitres1', fitres1
noise2 = fltarr(n_elements(x2)) + 1. ;set all noises equal
fitres2= MPFITFUN('linear',x2,y2, noise2, start)    ;ICL
;ap = plot(x1, y1)
;ap = plot(x1, fitres1(0)*x1 + fitres1(1),color = 'green',/overplot)
;a = plot(x2, y2, /overplot)
;ap = plot(x2, fitres2(0)*x2 + fitres2(1),color = 'green',/overplot)


flaty1 = y1 - (fitres1(0)*x1 + fitres1(1))
flaty2 = y2 - (fitres2(0)*x2+ fitres2(1))
;now put them together after flattening
x =[x1,x2]
flaty = [flaty1, flaty2]
;ap = plot(x, flaty+0.175, color = 'blue')

;now put them together after flattening
;x =[x1, x2]
;flaty = [flaty1, flaty2]
;print, y[0:100]
; Test the hypothesis that X and Y represent a significant periodic
; signal against the hypothesis that they represent random noise:
result = LNP_TEST(x, flaty,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
;PRINT, 'lnp result', result

b3 = plot(wk1, wk2, thick = 2, name = 'Self_cal', color = 'blue',/overplot);
;l = legend(target = [b2, b, b3, poly], position = [35, 140],/data)

end


