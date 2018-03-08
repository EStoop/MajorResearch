pro test_back_explosion_v2
;  ra_ref = 133.14754
;  dec_ref = 28.330195
  fwhm = 1.415

  ;ch2
  aorname = ['pcrs_planets/hd7924/r44605184', 'pcrs_planets/55cnc/r43981056', 'pcrs_planets/55cnc/r43981312', 'pcrs_planets/55cnc/r39524608', 'pcrs_planets/55cnc/r43981568', 'pcrs_planets/55cnc/r43981824', 'staring/r34769408', 'staring/r34869504', 'staring/r42000384', 'staring/r42790912', 'staring/r44273920']
  binsize = [64, 64, 100*64L,64L, 64, 64, 7, 7]

  ;ch1
;  aorname = ['pcrs_planets/wasp33/r39445760', 'pcrs_planets/wasp33/r45383680', 'pcrs_planets/wasp33/r45383936', 'pcrs_planets/wasp33/r45384192','staring/r34776832', 'staring/r38110464', 'staring/r34769152', 'staring/r38110720']
;  binsize = [7,7,7,7,7,7,64, 64]

  for a = 0, 0 do begin;  n_elements(aorname) - 1 do begin
     print, 'working on aorname ', aorname(a)
     dir = '/Users/jkrick/irac_warm/' + aorname(a) + '/ch2/bcd/'
     CD, dir                    ; change directories to the corrct AOR directory
     command  =  "find . -name '*_bcd.fits' > /Users/jkrick/irac_warm/staring/bcdlist.txt"
     spawn, command
     readcol,'/Users/jkrick/irac_warm/staring/bcdlist.txt',fitsname, format = 'A', /silent

     xarr = dblarr(64*n_elements(fitsname))
     yarr = dblarr(64*n_elements(fitsname))
     skyarr = dblarr(64*n_elements(fitsname))
     fluxarr = dblarr(64*n_elements(fitsname))
     backarr = dblarr(64*n_elements(fitsname))
     backsigarr = dblarr(64*n_elements(fitsname))
     sclkarr = dblarr(64*n_elements(fitsname))
     prfarr = dblarr(64*n_elements(fitsname))

     i = 0L
     print, 'n fits', n_elements(fitsname)
     ap3 = [3.]
     back3 = [11.,15.]           ; [11., 15.5]
     
     for f =0.D,  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        fits_read, fitsname(f), data, header  
        gain = sxpar(header, 'GAIN')
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        ronoise = sxpar(header, 'RONOISE')
        aorlabel = sxpar(header, 'OBJECT')
        sclk = sxpar(header, 'SCLK_OBS')

        convfac = gain*exptime/fluxconv
;        adxy, header, ra_ref, dec_ref, xinit, yinit
        xinit = 15.5
        yinit = 15.5
        for j = 0, 63 do begin
           indim = data[*,*,j]
           indim = indim*convfac
           irac_box_centroider, indim, xinit, yinit, 3, 6, 3, ronoise,xcen, ycen, box_f, box_sky, box_c, box_cb, box_np,/MMM
           xarr[i] = xcen
           yarr[i] = ycen
           skyarr[i] = box_sky
           
;now run aper at the centroids found, and use that background estimate.
           aper, indim,xcen, ycen, xf, xfs, xb, xbs, 1.0, ap3, back3, $
                 /FLUX, /EXACT, /SILENT, /NAN, READNOISE=ronoise
           fluxarr[i] = xf
           backarr[i] = xb
           backsigarr[i] =xbs
           sclkarr[i] = sclk

;try clacluating prf width
           aper, indim, xcen, ycen, topflux, topfluxerr, xb, xbs, 1.0, 8,[10,12],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
           aper, indim^2, xcen, ycen, bottomflux, bottomfluxerr, xb, xbs, 1.0,8,[10,12],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0

           beta = topflux^2 / bottomflux
           prfarr[i] = beta
;p----------------

           i = i + 1

        endfor
     endfor
     print, 'i', i
     xarr = xarr[0:i-1]
     yarr = yarr[0:i-1]
     skyarr = skyarr[0:i-1]
     fluxarr = fluxarr[0:i-1]
     backarr = backarr[0:i-1]
     backsigarr = backsigarr[0:i-1]
     sclkarr = sclkarr[0:i-1]
     numberarr = findgen(n_elements(xarr))
     prfarr = prfarr[0:i-1]
    
;  a = plot(numberarr, xarr,'1.', xtitle = 'number', ytitle = 'x cntrd', yrange = [14.4, 14.9])
;  a = plot(numberarr, yarr,'1.', xtitle = 'number', ytitle = 'y cntrd', yrange = [14.9, 15.1])
;     a = plot(numberarr, skyarr, '1.', xtitle = 'number', ytitle = 'box_sky')
;  a = plot(numberarr, fluxarr,  '1.', xtitle = 'number', ytitle = 'aper_flux')
;     a = plot(numberarr, backarr,  '1.', xtitle = 'number', ytitle = 'aper_bkgd')
     
;----------------------------
;try binning
;histogram numberarr by 64 as a technique to bin
     h = histogram(numberarr, OMIN=om, binsize = binsize(a), reverse_indices = ri)
                                ;print, 'h', h
     print, 'omin', om, 'nh', n_elements(h)
     
;mean together the flux values in each phase bin
     binxcen = dblarr(n_elements(h))
     binnum = dblarr(n_elements(h))
     binycen= dblarr(n_elements(h))
     binbkgd =  dblarr(n_elements(h))
     binflux=  dblarr(n_elements(h))
     binback=  dblarr(n_elements(h))
     bintime = dblarr(n_elements(h))
     binbacksig = dblarr(n_elements(h))
     binprf = dblarr(n_elements(h))

     c = 0L
     for j = 0L, n_elements(h) - 2 do begin
        
;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j])  then begin
                                ;print, 'binning together', n_elements(numberarr[ri[ri[j]:ri[j+1]-1]])
                                ;print, 'binning', numberarr[ri[ri[j]:ri[j+1]-1]]
           
           meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           binxcen[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
           binycen[c] = meany   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
;           meanclip, skyarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
;           binbkgd[c] = meansky ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           binnum[c] = mean(numberarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, fluxarr[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           binflux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, backarr[ri[ri[j]:ri[j+1]-1]], meanback, sigmaback
           mmm, backarr[ri[ri[j]:ri[j+1]-1]],skymod, skysigma, skew, /silent
           ;print, 'meanclip, mmm', sigmaback, skysigma
           binback[c] = skymod; meanback ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, backsigarr[ri[ri[j]:ri[j+1]-1]], meanbacksig, sigmabacksig

           binbacksig[c]= skysigma; sigmaback

           bintime[c] =mean(sclkarr[ri[ri[j]:ri[j+1]-1]])

           ;plot some histograms
          ; plothist,backarr[ri[ri[j]:ri[j+1]-1]], xhist, yhist, bin = 0.4, /noplot
          ; if c eq 0 then begin
          ;    sp = plot(xhist, yhist, xtitle = 'background (electrons)') 
          ; endif else begin
          ;    sp = plot(xhist, yhist - j*50., /overplot)
         ;  endelse
          
           meanclip, prfarr[ri[ri[j]:ri[j+1]-1]], meanprf, sigmaprf
           binprf[c] = meanprf 

           c = c + 1
                                ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
        endif
     endfor
     
     binxcen = binxcen[0:c-1]
     binycen = binycen[0:c-1]
     binnum = binnum[0:c-1]
     binbkgd = binbkgd[0:c-1]
     bintime = bintime[0:c-1]
     binflux = binflux[0:c-1]
     binback = binback[0:c-1]
     binbacksig = binbacksig[0:c-1]
     binprf = binprf[0:c-1]

     save, /all, filename = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/back_save_bin.sav'
     
     print, 'c', c
                                ;and plot
;  a = plot(binnum, binxcen, '1.', xtitle = 'number', ytitle = 'binned x cntrd', yrange = [14.4, 14.9])
;  a = plot(binnum, binycen, '1.', xtitle = 'number', ytitle = 'binned y cntrd')
;     a = plot(binnum, binbkgd, '1.', xtitle = 'number', ytitle = 'binned box_bkgdd')
;  a = plot(binnum, binflux, '1.', xtitle = 'number', ytitle = 'binned aper_flux')
     titlename = string(aorname(a)) +' '+ string(aorlabel) + string(exptime)
     b = plot((bintime - bintime(0))/60./60., binback, '1', xtitle = 'Time (hrs)', ytitle = 'binned background (electrons)', title = titlename)

     


;now for the correlation plots:
p = plot((bintime - bintime(0)) / 60./60., binbacksig, '1s', sym_size = 0.1, ytitle = 'background sigma',xtitle = 'time(hrs)')

p = plot(binback, binbacksig, '1s',sym_size = 0.1, xtitle = 'binned background (electrons)', ytitle = 'background sigma')

p = plot((bintime - bintime(0)) / 60./60., binprf, '1s', sym_size = 0.1, ytitle = 'prf size',xtitle = 'time(hrs)')
p = plot((bintime - bintime(0)) / 60./60., binxcen, '1s', sym_size = 0.1, ytitle = 'xcen',xtitle = 'time(hrs)')
p = plot((bintime - bintime(0)) / 60./60., binycen, '1s', sym_size = 0.1, ytitle = 'ycen',xtitle = 'time(hrs)')

endfor
end
 
