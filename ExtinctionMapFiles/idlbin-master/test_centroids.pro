pro test_centroids
  ra_ref = 133.14754
  dec_ref = 28.330195
  fwhm = 1.415
  
  dir = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/r39524608/ch2/bcd/'
  CD, dir                       ; change directories to the corrct AOR directory
  command  =  "find . -name '*_bcd.fits' > /Users/jkrick/irac_warm/pcrs_planets/55cnc/bcdlist.txt"
  spawn, command
  readcol,'/Users/jkrick/irac_warm/pcrs_planets/55cnc/bcdlist.txt',fitsname, format = 'A', /silent

  xarr = dblarr(64*n_elements(fitsname))
  yarr = dblarr(64*n_elements(fitsname))
  skyarr = dblarr(64*n_elements(fitsname))
  fluxarr = dblarr(64*n_elements(fitsname))
  backarr = dblarr(64*n_elements(fitsname))
  i = 0L
  print, 'n fits', n_elements(fitsname)
  ap3 = [3.]
  back3 = [11.,15.]; [11., 15.5]

  for f =0.D,  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
     fits_read, fitsname(f), data, header  
     gain = sxpar(header, 'GAIN')
     fluxconv = sxpar(header, 'FLUXCONV')
     exptime = sxpar(header, 'EXPTIME')
     ronoise = sxpar(header, 'RONOISE')
     convfac = gain*exptime/fluxconv
     adxy, header, ra_ref, dec_ref, xinit, yinit
     for j = 0, 63 do begin
        indim = data[*,*,j]
        indim = indim*convfac
        if f eq 10 and j eq 1 then fits_write, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/r39524608/ch2/bcd/electrons.fits', indim, header
        irac_box_centroider, indim, xinit, yinit, 3, 6, 3, ronoise,xcen, ycen, box_f, box_sky, box_c, box_cb, box_np,/MMM
        xarr[i] = xcen
        yarr[i] = ycen
        skyarr[i] = box_sky
        
;now run aper at the centroids found, and use that background estimate.
        aper, indim,xcen, ycen, xf, xfs, xb, xbs, 1.0, ap3, back3, $
              /FLUX, /EXACT, /SILENT, /NAN, READNOISE=ronoise
              
        fluxarr[i] = xf
        backarr[i] = xb
        
        i = i + 1
     endfor
  endfor
  print, 'i', i
  xarr = xarr[0:i-1]
  yarr = yarr[0:i-1]
  skyarr = skyarr[0:i-1]
  fluxarr = fluxarr[0:i-1]
  backarr = backarr[0:i-1]
  numberarr = findgen(n_elements(xarr))
  
  
  
;  a = plot(numberarr, xarr,'1.', xtitle = 'number', ytitle = 'x cntrd', yrange = [14.4, 14.9])
;  a = plot(numberarr, yarr,'1.', xtitle = 'number', ytitle = 'y cntrd', yrange = [14.9, 15.1])
;  a = plot(numberarr, skyarr, '1.', xtitle = 'number', ytitle = 'box_sky')
;  a = plot(numberarr, fluxarr,  '1.', xtitle = 'number', ytitle = 'aper_flux')
;  a = plot(numberarr, backarr,  '1.', xtitle = 'number', ytitle = 'aper_bkgd')

;----------------------------
;try binning
;histogram numberarr by 64 as a technique to bin
    h = histogram(numberarr, OMIN=om, binsize = 64, reverse_indices = ri)
    ;print, 'h', h
    print, 'omin', om, 'nh', n_elements(h)

;mean together the flux values in each phase bin
  binxcen = dblarr(n_elements(h))
  binnum = dblarr(n_elements(h))
  binycen= dblarr(n_elements(h))
  binbkgd =  dblarr(n_elements(h))
  binflux=  dblarr(n_elements(h))
  binback=  dblarr(n_elements(h))
  c = 0
  for j = 0L, n_elements(h) - 1 do begin

;get rid of the bins with no values and low numbers, meaning low overlap
     if (ri[j+1] gt ri[j])  then begin
        ;print, 'binning together', n_elements(numberarr[ri[ri[j]:ri[j+1]-1]])
        ;print, 'binning', numberarr[ri[ri[j]:ri[j+1]-1]]
        
        meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        binxcen[c] = meanx; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
        binycen[c] = meany; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, skyarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
        binbkgd[c] = meansky; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        binnum[c] = mean(numberarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, fluxarr[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
        binflux[c] = meanflux; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, backarr[ri[ri[j]:ri[j+1]-1]], meanback, sigmaback
        binback[c] = meanback; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        c = c + 1
        ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
     endif
  endfor

  binxcen = binxcen[0:c-1]
  binycen = binycen[0:c-1]
  binnum = binnum[0:c-1]
  binbkgd = binbkgd[0:c-1]
;  binflux = binflux[0:c-1]
;  binback = binback[0:c-1]
  
  print, 'c', c
                                ;and plot
;  a = plot(binnum, binxcen, '1.', xtitle = 'number', ytitle = 'binned x cntrd', yrange = [14.4, 14.9])
;  a = plot(binnum, binycen, '1.', xtitle = 'number', ytitle = 'binned y cntrd')
;  a = plot(binnum, binbkgd, '1.', xtitle = 'number', ytitle = 'binned box_bkgdd')
;  a = plot(binnum, binflux, '1.', xtitle = 'number', ytitle = 'binned aper_flux')
  a = plot(binnum, binback, '1.', xtitle = 'number', ytitle = 'binned aper_back')

end
