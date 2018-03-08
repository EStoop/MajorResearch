pro plot_wasp33_ch2
planetname = 'WASP-33b'
  colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']

;  fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp-33b/WASP-33b_phot_ch2_2.25000superdark.sav'
  aorname = ['r45383424', 'r45384448', 'r45384704'] ;ch2
  chname = '2'

;-------------------------------------------------------
;need to re-correct with the new pmap
; AORwasp33 = replicate({wasp33ob, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits)},n_elements(aorname))

;  file_suffix = ['500x500_0043_120409.fits','0p1s_x4_500x500_0043_120531.fits']
;print, 'starting corrflux ', systime()
;print, AORwasp33[0].corrflux[10:100]
;  AORwasp33.corrflux = iracpc_pmap_corr(AORwasp33.flux,AORwasp33.xcen,AORwasp33.ycen,2,FILE_SUFFIX=file_suffix,/threshold_occ, threshold_val = 20)
;print, 'ending corrflux ', systime()
;print, AORwasp33[0].corrflux[10:100]

;  save, AORwasp33, filename='/Users/jkrick/irac_warm/pcrs_planets/wasp33/wasp33_phot_ch2_pmap2.sav'


;-------------------------------------------------------
;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
;  chname = planetinfo[planetname, 'chname']
;  ra_ref = planetinfo[planetname, 'ra']
;  dec_ref = planetinfo[planetname, 'dec']
  aorname = planetinfo[planetname, 'aorname']
  basedir = planetinfo[planetname, 'basedir']
;  utmjd_center =  planetinfo[planetname, 'utmjd_center']
;  transit_duration =  planetinfo[planetname, 'transit_duration']
;  period =  planetinfo[planetname, 'period']
;  intended_phase = planetinfo[planetname, 'intended_phase']
;  stareaor = planetinfo[planetname, 'stareaor']
;  plot_norm= planetinfo[planetname, 'plot_norm']
;  plot_corrnorm = planetinfo[planetname, 'plot_corrnorm']
 


;
; c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500])

;  for a = 0, n_elements(aorname) - 1 do begin
;     xcen500 = 500.* (AORwasp33[a].xcen - 14.5)
;     ycen500 = 500.* (AORwasp33[a].ycen - 14.5)
;     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
; endfor

;-----
 
;    for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'X pix', yrange = [14.8, 15.3], xrange = [0,40])
;     endif else begin
;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------

;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Y pix',  yrange = [14.8, 15.3], xrange = [0,40])
;     endif else begin
;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Background', xrange = [0,40])
;     endif else begin
;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse
;
;  endfor

;------

 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;try this diffferently, try putting the AORs togetether as one array; then bin

timearr = [planethash[aorname(0),'timearr'], planethash[aorname(1),'timearr'], planethash[aorname(2),'timearr']]
fluxarr = [planethash[aorname(0),'flux'], planethash[aorname(1),'flux'], planethash[aorname(2),'flux']]
corrfluxarr = [planethash[aorname(0),'corrflux'], planethash[aorname(1),'corrflux'], planethash[aorname(2),'corrflux']]
corrfluxerrarr =  [planethash[aorname(0),'corrfluxerr'], planethash[aorname(1),'corrfluxerr'], planethash[aorname(2),'corrfluxerr']]
xarr = [planethash[aorname(0),'xcen'], planethash[aorname(1),'xcen'], planethash[aorname(2),'xcen']]
yarr = [planethash[aorname(0),'ycen'], planethash[aorname(1),'ycen'], planethash[aorname(2),'ycen']]
bkgd = [planethash[aorname(0),'bkgd'], planethash[aorname(1),'bkgd'], planethash[aorname(2),'bkgd']]
bmjd = [planethash[aorname(0),'bmjdarr'] , planethash[aorname(1),'bmjdarr'] , planethash[aorname(2),'bmjdarr'] ]

print, 'n', n_elements(timearr), n_elements(planethash[aorname(0),'timearr']), n_elements(planethash[aorname(1),'timearr']), n_elements(planethash[aorname(2),'timearr'])

timearr = timearr[0:103005];[206010:*];[103005:206010];
fluxarr = fluxarr[0:103005];[206010:*];[103005:206010];
corrfluxarr = corrfluxarr[0:103005];[206010:*];[103005:206010];
xarr = xarr[0:103005];[206010:*];[103005:206010];
yarr = yarr[0:103005];[206010:*];[103005:206010];[
bkgd = bkgd[0:103005];[206010:*];[103005:206010];[
bmjd = bmjd[0:103005];[206010:*];[103005:206010];


;now for the binning

      ;remove outliers
good = where(xarr lt mean(xarr) + 2.5*stddev(xarr) and xarr gt mean(xarr) -2.5*stddev(xarr) and xarr lt mean(xarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr),ngood_pmap, complement=bad) 
     print, 'bad ',n_elements(bad)
     xarr = xarr[good]
     yarr = yarr[good]
     timearr = timearr[good]
     flux = fluxarr[good]
     corrflux = corrfluxarr[good]
     corrfluxerr = corrfluxerrarr[good]
     bmjdarr = bmjd[good]
     bkgdarr = bkgd[good]
     
     
; binning
     bin_level = 63L;*60
     nframes = bin_level
     nfits = long(n_elements(corrflux)) / nframes
     print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
     bin_flux = dblarr(nfits)
     bin_corrflux= dblarr(nfits)
     bin_corrfluxerr= dblarr(nfits)
     bin_timearr = dblarr(nfits)
     bin_bmjdarr = dblarr(nfits)
     bin_bkgd = dblarr(nfits)

     for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', n_elements(AORwasp33[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
        idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
        id = flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
        ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
       ; if si lt 2 then print, 'ib', ib
       ; print, 'mean', mean(ib,/nan)
        bin_flux[si] = mean(id, /nan)
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
     bin_corrflux = bin_corrflux[good_ni]
     bin_bkgd = bin_bkgd[good_ni]

     print, 'nflux', n_elements(bin_flux)
     
;now try plotting
save, bin_timearr, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_timearr.sav'
save, bin_flux, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_flux.sav'
save, bin_corrflux, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_corrflux.sav'

     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux, '6r1s', sym_size = 0.1,   sym_filled = 1, xrange = [0,40], yrange = [0.165, 0.180], xtitle = 'Time (hrs)', ytitle = 'Flux (Jy)')
     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux-0.002, '6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)

;     pl = plot(bin_bmjdarr, bin_flux/median(bin_flux) + 0.03, '6r1s', sym_size = 0.1,   sym_filled = 1,  yrange = [0.97,1.05], xrange = [-0.7,0.7],xtitle = 'Phase', ytitle = 'Normalized Flux')
;     pl = plot(bin_bmjdarr, (bin_corrflux/median(bin_corrflux))-0.003, '6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)

;------------------------------------------------------------
;------------------------------------------------------------
;periodogram
;------------------------------------------------------------
;------------------------------------------------------------
;print, 'starting periodogram'
;X = (bin_timearr - bin_timearr(0))/60./60. ;in hours
X = (bin_timearr - bin_timearr(0))/60./60./24. ;in days

Y = bin_corrflux
;t = plot(x, y)
;ok, try cutting this down in phase to between eclipses
;a = where(x gt 7 and x lt 17)  ;in hours
a1 = where(x gt 0.25 and x lt 0.7) ; in days
a2 = where(x gt 1.08 and x lt 1.33) ; in days 0.9

x1 =x(a1)
x2 = x(a2)
y1 = y(a1)
y2=y(a2)
;badx = where(finite(x) lt 1, badxcount)
bady = where(finite(y1) lt 1, badycount)
print,'badycount', badycount
;print, x(bady), y(bady)
if badycount gt 0 then y1(bady) = y1(bady-1) 
;print, x(bady), y(bady)
;print, y[0:100]
ytot = [y1, y2]

;set the mean to zero
;doesn't seem to make a difference.
;y = y - mean(y)

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
;
flaty1 = y1 - (fitres1(0)*x1 + fitres1(1))
flaty2 = y2 - (fitres2(0)*x2+ fitres2(1))
;ap = plot(x, flaty, color = 'blue', /overplot)

;now put them together after flattening
x = [x1, x2]
flaty = [flaty1, flaty2]

;convert x to minutes
;x = x / 24. / 60.

print, 'pmap units', min(x), max(x), min(flaty), max(flaty)
;print, y[0:100]
; Test the hypothesis that X and Y represent a significant periodic
; signal against the hypothesis that they represent random noise:
result = LNP_TEST(X, flaty,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
PRINT, 'lnp result', result

b = plot(wk1, wk2, xtitle = 'Frequency(1/days)', ytitle = 'Power', xrange =[0,50], yrange = [0,150],thick = 2, color = 'red',name = 'Pmap Corrected')

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
print, 'uncorrected units', min(x), max(x), min(flaty2), max(flaty2)

result = LNP_TEST(X, flatY2,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
PRINT, result

b2 = plot(wk1, wk2, /overplot, color = 'black',thick = 2, name = 'Original')
save, wk1, wk2, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp33/periodogram1.sav'
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

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_time.sav'
restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_sub.sav'
x = bin_time /24.
y = bin_sub +0.175

a1 = where(x gt 6.5/24. and x lt 17.5/24.)
x1 =x(a1)
print, 'fraction using' , n_elements(x1), n_elements(x)
y1 = y(a1)
a2 = where(x gt 21.1/24. and x lt 31.5/24.)
print, 'fraction using' , n_elements(x2), n_elements(x)
x2 = x(a2)
y2 = y(a2)

;badx = where(finite(x1) lt 1, badxcount)
;bady = where(finite(y2) lt 1, badycount)
;print,'badycount', badycount, badxcount
;if badycount gt 0 then y2(bady) = y2(bady-1) 


;try taking out the phase variation
;which also takes out the mean
start =  [-0.001,0.173]
noise1 = fltarr(n_elements(x1)) + 1. ;set all noises equal
fitres1= MPFITFUN('linear',x1,y1, noise1, start)    ;ICL
print, 'fitres1', fitres1
;ap = plot(x1, y1)
;ap = plot(x1, fitres1(0)*x1 + fitres1(1),color = 'green',/overplot)
noise2 = fltarr(n_elements(x2)) + 1. ;set all noises equal
fitres2= MPFITFUN('linear',x2,y2, noise2, start)    ;ICL
print, 'fitres1', fitres1

flaty1 = y1 - (fitres1(0)*x1 + fitres1(1))
flaty2 = y2 - (fitres2(0)*x2 + fitres2(1))
x = [x1,x2]
flaty = [flaty1,flaty2]
;ap = plot(x, flaty+0.175, color = 'blue')

;now put them together after flattening
; Test the hypothesis that X and Y represent a significant periodic
; signal against the hypothesis that they represent random noise:
print, 'selfcal units', min(x), max(x), min(flaty), max(flaty)

result = LNP_TEST(x, flaty,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
PRINT, 'lnp result', result

;b3 = plot(wk1, wk2, thick = 2, name = 'Self-Calibrated', color = 'light_sky_blue',/overplot);
;l = legend(target = [b2, b, b3], position = [35, 140],/data)

end


;old stuff

; for a = 0, n_elements(aorname) - 1 do begin
;;sort by time
;     s = sort(AORwasp33[a].bmjdarr)
;     sorttime = AORwasp33[a].timearr[s]
;     sortxcen= AORwasp33[a].xcen[s]
;     sortycen = AORwasp33[a].ycen[s]
;     sortflux = AORwasp33[a].flux[s]
;     sortcorrflux = AORwasp33[a].corrflux[s]
;     sortcorrfluxerr = AORwasp33[a].corrfluxerr[s]
;     sortbmjd = AORwasp33[a].bmjdarr[s]
;     sortutcs = AORwasp33[a].utcsarr[s]
;     sortbkgd = AORwasp33[a].bkgd[s]
;
;      ;remove outliers
;     good = where(sortxcen lt mean(sortxcen) + 2.5*stddev(sortxcen) and sortxcen gt mean(sortxcen) -2.5*stddev(sortxcen) and sortycen lt mean(sortycen) +3.0*stddev(sortycen) and sortycen gt mean(sortycen) - 3.0*stddev(sortycen),ngood_pmap, complement=bad) 
;     print, 'bad ',n_elements(bad)
;     xarr = sortxcen[good]
;     yarr = sortycen[good]
;     timearr = sorttime[good]
;     flux = sortflux[good]
;     corrflux = sortcorrflux[good]
;     corrfluxerr = sortcorrfluxerr[good]
;     bmjdarr = sortbmjd[good]
;     bkgdarr = sortbkgd[good]
;     
;     
;; binning
;     bin_level = 63L;*60
;     nframes = bin_level
;     nfits = long(n_elements(corrflux)) / nframes
;     print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
;     bin_flux = dblarr(nfits)
;     bin_corrflux= dblarr(nfits)
;     bin_corrfluxerr= dblarr(nfits)
;     bin_timearr = dblarr(nfits)
;     bin_bmjdarr = dblarr(nfits)
;     bin_bkgd = dblarr(nfits)
;
;     for si = 0L, long(nfits) - 1L do begin
;                                ;print, 'working on si', n_elements(AORwasp33[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
;        idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
;        id = flux[si*nframes:si*nframes + (nframes - 1)]
;        idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
;        ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
;       ; if si lt 2 then print, 'ib', ib
;       ; print, 'mean', mean(ib,/nan)
;        bin_flux[si] = mean(id, /nan)
;        bin_corrflux[si] = mean(idata,/nan)
;        bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
;        bin_timearr[si]= mean(timearr[si*nframes:si*nframes + (nframes - 1)])
;        bin_bmjdarr[si]= mean(bmjdarr[si*nframes:si*nframes + (nframes - 1)])
;        bin_bkgd[si] = mean(ib,/nan)
;
;     endfor                     ;for each fits image
;     
; ;try getting rid of flux outliers.
;  ;do some running mean with clipping
;     start = 0
;     print, 'nflux', n_elements(bin_flux)
;     for ni = 50, n_elements(bin_flux) -1, 50 do begin
;        meanclip,flux[start:ni], m, s, subs = subs ;,/verbose
;                                ;print, 'good', subs+start
;                                ;now keep a list of the good ones
;        if ni eq 50 then good_ni = subs+start else good_ni = [good_ni, subs+start]
;        start = ni
;     endfor
;                                ;see if it worked
;     ;xarr = xarr[good_ni]
;     ;yarr = yarr[good_ni]
;     bin_timearr = bin_timearr[good_ni]
;     bin_bmjdarr = bin_bmjdarr[good_ni]
;     bin_corrfluxerr = bin_corrfluxerr[good_ni]
;     bin_flux = bin_flux[good_ni]
;     bin_corrflux = bin_corrflux[good_ni]
;     bin_bkgd = bin_bkgd[good_ni]
;
;     print, 'nflux', n_elements(bin_flux)
;     
; 
;;     if a eq 0 then begin
;;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Flux', xrange = [0,40], yrange = [0.16,0.18])
;;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].corrflux-0.01,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)
;
;;     endif else begin
;;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;;        am = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].corrflux-0.01,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)
;;     endelse
;
;
;  endfor
;  


; and the periodogram for the self calibrated data

;restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_time.sav'
;restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_sub.sav'
;x = bin_time 
;y = bin_sub
;
;a1 = where(x gt 0 and x lt 3.0) ; in days
;a2 = where(x gt 6.5 and x lt 12.3)
;a3 = where(x gt 12.4 and x lt 17.5) ; in days 0.9
;a4 = where(x gt 21. and x lt 24.7)
;a5 = where(x gt 25 and x lt 32.)
;
;a1 = where(x gt 6.5 and x lt 17.5)
;
;x1 =x(a1)
;print, 'fraction using' , n_elements(x1), n_elements(x)
;x2 = x(a2)
;x3 = x(a3)
;x4 = x(a4)
;x5 = x(a5)
;
;y1 = y(a1)
;y2=y(a2)
;y3=y(a3)
;y4=y(a4)
;y5=y(a5)
;
;;badx = where(finite(x) lt 1, badxcount)
;;bady = where(finite(y1) lt 1, badycount)
;;print,'badycount', badycount
;;;print, x(bady), y(bady)
;;if badycount gt 0 then y1(bady) = y1(bady-1) 
;;;print, x(bady), y(bady)
;;;print, y[0:100]
;;ytot = [y1, y2,y3,y4,y5]
;
;;try taking out the phase variation
;;which also takes out the mean
;start =  [-0.001,0.173]
;noise1 = fltarr(n_elements(x1)) + 1. ;set all noises equal
;fitres1= MPFITFUN('linear',x1,y1, noise1, start)    ;ICL
;print, 'fitres1', fitres1
;noise2 = fltarr(n_elements(x2)) + 1. ;set all noises equal
;fitres2= MPFITFUN('linear',x2,y2, noise2, start)    ;ICL
;noise3 = fltarr(n_elements(x3)) + 1. ;set all noises equal
;fitres3= MPFITFUN('linear',x3,y3, noise3, start)    ;ICL
;noise4 = fltarr(n_elements(x4)) + 1. ;set all noises equal
;fitres4= MPFITFUN('linear',x4,y4, noise4, start)    ;ICL
;noise5 = fltarr(n_elements(x5)) + 1. ;set all noises equal
;fitres5= MPFITFUN('linear',x5,y5, noise5, start)    ;ICL
;;ap = plot(x1, y1)
;;ap = plot(x1, fitres1(0)*x1 + fitres1(1),color = 'green',/overplot)
;;a = plot(x2, y2, /overplot)
;;ap = plot(x2, fitres2(0)*x2 + fitres2(1),color = 'green',/overplot)
;;a = plot(x3, y3, /overplot)
;;ap = plot(x3, fitres3(0)*x3 + fitres3(1),color = 'green',/overplot)
;;a = plot(x4, y4, /overplot)
;;ap = plot(x4, fitres4(0)*x4 + fitres4(1),color = 'green',/overplot)
;;a = plot(x5, y5, /overplot)
;;ap = plot(x5, fitres5(0)*x5 + fitres5(1),color = 'green',/overplot)
;
;
;flaty1 = y1 - (fitres1(0)*x1 + fitres1(1))
;flaty2 = y2 - (fitres2(0)*x2+ fitres2(1))
;flaty3 = y3 - (fitres3(0)*x3+ fitres3(1))
;flaty4 = y4 - (fitres4(0)*x4+ fitres4(1))
;flaty5 = y5 - (fitres5(0)*x5+ fitres5(1))
;;now put them together after flattening
;x =[x1,x2,x3,x4,x5]
;flaty = [ flaty1, flaty2,flaty3, flaty4, flaty5]
;
;x = x1
;flaty = flaty1
;;ap = plot(x, flaty+0.175, color = 'blue')
;
;;now put them together after flattening
;;x =[x1, x2]
;;flaty = [flaty1, flaty2]
;;print, y[0:100]
;; Test the hypothesis that X and Y represent a significant periodic
;; signal against the hypothesis that they represent random noise:
;result = LNP_TEST(x, flaty,/double, WK1 = wk1, WK2 = wk2, JMAX = jmax)
;;PRINT, 'lnp result', result
;
;b3 = plot(wk1, wk2, thick = 2, name = 'Self_cal', color = 'red',/overplot);
;;l = legend(target = [b2, b, b3, poly], position = [35, 140],/data)
;
