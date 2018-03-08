function color_by_exptime, exptime

  case exptime of
     0.01: plotcolor ='red'
     0.08: plotcolor = 'coral'
     0.36: plotcolor = 'dark_orange'
     1.92: plotcolor = 'yellow'
     4.40: plotcolor = 'chartreuse'
     10.4: plotcolor = 'dark_green'
     26.8: plotcolor = 'cyan'
     96.8: plotcolor =  'navy'
     else: plotcolor = 'black'
  endcase

  return, plotcolor
end


function find_binscale, exptime
  ;;make the largest binscale 1/3 the number of frames
  case exptime of
    
     0.01: binscale = fix(90.*60./0.01 /3.)
     0.08: binscale = fix(90.*60./0.1 /3.)
     0.36: binscale = fix(90.*60./0.4 /3.)
     1.92: binscale = fix(90.*60./2.0 /3.)
     4.40: binscale = fix(90.*60./6.0 /3.)
     10.4: binscale = fix(90.*60./12.0 /3.)
     26.8: binscale = fix(90.*60./30.0 /3.)
     96.8: binscale = fix(90.*60./100.0 /3.)
     else: binscale = fix(90.*60./2.0 /3.)
  endcase

  return, binscale
end

pro rms_binsize_alldata
  
  savenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08.sav']
  
  extra={xthick: 2, ythick:2, margin: 0.2, sym_size: 0.2,   sym_filled: 1, xstyle:1}
  
  starts = 0
  stops =  7                   ; n_elements(savenames) - 1
  totalaorcount = 0
  TIC

  for s = starts, stops do begin
     print, 'restoring ', savenames(s)
     restore, savenames(s)
     aorname = planethash.keys()
;;     allrms =fltarr(n_elements(aorname), bin_scale )
     
     for a = 0,  n_elements(aorname) - 1 do begin
        print, '--------------------------------'
        print, 'working on ', aorname(a), ' ', a
        corrflux =planethash[aorname(a),'pld_flux'] ;'corrflux_d'
        ;;corrfluxerr =planethash[aorname(a),'pld_unc'] ;'corrfluxerr'
;;        corrfluxerr = corrflux - corrflux + 1.0  ;fake it for now
        ;;phasearr = planethash[aorname(a),'phase']
        bmjdarr = planethash[aorname(a), 'bmjdarr']
        bmjdarr = bmjdarr - bmjdarr(0)
        xcen = planethash[aorname(a), 'xcen']
        ycen = planethash[aorname(a), 'ycen']
        min_dur = planethash[aorname(a), 'min_dur']
        exptime = planethash[aorname(a),'exptime']
        if exptime eq 0.01 then exptime = 0.05 ; otherwise total times of observation don't work out
        chname = planethash[aorname(a),'chname']
        flux = planethash[aorname(a), 'flux']
        ;; print, 'stdddev, min_dur', stddev(xcen), stddev(ycen), min_dur, exptime, chname
        
        ;;remove dithered observations and short AORs
        if (stddev(xcen) lt 0.15) and (stddev(ycen) lt 0.15) and (min_dur(n_elements(min_dur) - 1) gt 130) then begin
           ;;print, 'made the cut'
           
           ;;want to cut out the first 30 minutes in case there is
           ;;drift and only keep the first 2hrs to possibly get the
           ;;cleanest part of the light curve with no transits/eclipses
           cutstart = long(30.*60./ exptime) ;30min*60s/exptime = number of photometry points to remove up front
           cutend = long(120.*60./exptime)
           print, 'cutstart, cutend', cutstart, cutend, exptime, n_elements(corrflux)
           ;;some min_durs are incorrect?
           if cutend lt n_elements(corrflux) then begin
              corrflux = corrflux[cutstart:cutend]
              flux = flux[cutstart:cutend]
              print, 'after cutting light curve', exptime,  n_elements(corrflux)
              
              ;;get rid of NANs
              ;; b = where(finite(corrflux) gt 0, good) ;and phasearr lt -0.002 and phasearr gt -0.012
              ;;corrflux = corrflux(b)
              ;;corrfluxerr = corrfluxerr(b)
              ;; phasearr = phasearr(b)
              ;;bmjdarr = bmjdarr(b)
              
              ;;figure out the bin_scale based on exptime
              bin_scale = find_binscale(exptime)
              print, 'bin_scale', bin_scale
              ;;set up arrays 
              rmsarr_raw = fltarr(bin_scale + 1 )
              numberarr = findgen(n_elements(flux)) ; should be the same as PLD
              rmsarr_pld = fltarr(bin_scale + 1 )
              
              ;;color code by exptime
              plotcolor = color_by_exptime(exptime)
              
              for n = 1, bin_scale - 1  do begin
                 
                 ;;bin up the data
                 h = histogram(numberarr, OMIN=om, binsize = n, reverse_indices = ri)
                 
                 ;;setup some arrays
                 bin_corrflux = dblarr(n_elements(h))
                 bin_timearr = bin_corrflux
                 bin_flux = bin_corrflux
                 
                 if n eq 1 then begin                       ; no binning
                    rmsarr_pld(n) = robust_sigma(bin_corrflux) ; stddev(corrflux,/nan)
                    rmsarr_raw(n) = robust_sigma(flux)
                 endif else begin ; for the binning scales
                    c = 0
                    for j = 0L, n_elements(h) - 1 do begin
                       
                       if (ri[j+1] gt ri[j] + 1)  then begin ;require 2 elements in the bin
                          ;;what is the mean flux in each bin
                          meanclip, corrflux[ri[ri[j]:ri[j+1]-1]], meancorrflux, sigmacorrflux
                          bin_corrflux[c] = meancorrflux
                          
                          meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
                          bin_flux[c] = meanflux
                          
                          meantimearr = median(bmjdarr[ri[ri[j]:ri[j+1]-1]])
                          bin_timearr[c]=meantimearr
                          
                          c = c + 1
                       endif
                    endfor
                    bin_corrflux = bin_corrflux[0:c-1]
                    bin_flux = bin_flux[0:c-1]
                    bin_timearr = bin_timearr[0:c-1]
                    
                    meanclip, bin_corrflux, robustmean, robustsigma
                    meanclip, bin_flux, robustmeanraw, robustsigmaraw
                    rmsarr_pld(n) =  robustsigma ; robust_sigma(bin_corrflux)
                    rmsarr_raw(n) =  robustsigmaraw ; robust_sigma(bin_flux)
                                ;rmsarr(n) =  stddev(bin_corrflux,/nan)
                    
                 endelse
;;              allrms(a,n) = rmsarr(n)
                 
              endfor            ; for all binning scales
              
              if chname eq 'ch2' then begin
                 ;;p2 = plot(findgen(n_elements(rmsarr_pld)), rmsarr_pld/ rmsarr_pld(2), /xlog, /ylog, overplot = p2, $ ;
              ;;          color = plotcolor,   xtitle = 'Number of frames',yrange = [1E-2,1.5], $
                 ;;          ytitle = 'Normalized Residual RMS', title = 'CH2') ; axis_style = 1 ;xrange =[1, bin_scale],
              endif
              if chname eq 'ch1' then begin
                 ;; p1 = plot(findgen(n_elements(rmsarr_pld)), rmsarr_pld/ rmsarr_pld(2), /xlog, /ylog, overplot = p1, $ ;
                 ;;           color = plotcolor,  xtitle = 'Number of frames',yrange = [1E-2,1.5], $
                 ;;           ytitle = 'Normalized Residual RMS', title = 'CH1') ; axis_style = 1 xrange =[1, bin_scale],
              endif
              
              ;;measure the slope of the rmsarr starting at 10th point
              rmsarr_pld = rmsarr_pld[10:n_elements(rmsarr_pld) - 2]
              x = findgen(n_elements(rmsarr_pld)) + 10 ; since we chopped 0 - 10
              logx = alog10(x)
              logrmsarr_pld = alog10(rmsarr_pld)
              
              result = ladfit(logx, logrmsarr_pld, absdev = absdev_pld)
              slope_pld = result(1)
              
              rmsarr_raw = rmsarr_raw[10:n_elements(rmsarr_raw) - 2]
              logrmsarr_raw = alog10(rmsarr_raw)
              
              result = ladfit(logx, logrmsarr_raw, absdev = absdev_raw)
              slope_raw = result(1)
              print, 'slope raw, pld', slope_raw,  slope_pld, absdev_raw, absdev_pld
              
              ;;add data to the planethashs
              planethash[aorname(a),'rms_slope_pld'] = slope_pld
              planethash[aorname(a),'rms_absdev_pld'] = absdev_pld
              planethash[aorname(a),'rms_slope_raw'] = slope_raw
              planethash[aorname(a),'rms_absdev_raw'] = absdev_raw
              
              
           endif
        endif
        
     endfor                     ; for all AORs
     save, planethash, filename=savenames(s)
     
  endfor                        ; for all save files
  
  ;;what is the average or median of the allrms array?
 ;; medianrms = median(allrms, dimension = 1)
  
  ;;try normalizing
  ;;medianrms = medianrms/medianrms(1)
  ;;p2 = plot(findgen(n_elements(rmsarr)), medianrms, thick =3,color = 'black', overplot = p2)
;------------------------
;now add straight root N
     
  ;;source_mjy = 68.29            ;mJy  from Star-pet
  ;;exptime = 2.
     
;ch2
  ;;gain = 3.71
  ;;pixel_scale = 1.22
  ;;flux_conv = .1469
  ;;xmax = max(p2.xrange)
  ;;bin_scale = findgen(xmax) + 1
  ;;root_n = (bin_scale)^(0.5)
  ;;source_electrons = mjy_to_electron( source_mjy, pixel_scale, gain, exptime, flux_conv)
  ;;sigma_poisson = sqrt(source_electrons)
  ;;y =  (sigma_poisson / root_n) / source_electrons
  
  ;;p2 = plot( bin_scale, y/y(0), thick = 3,  overplot = p2)  ;
  ;;p1 = plot( bin_scale, y/y(0), thick = 3,  overplot = p1)  ;
;now add root2 * root N
  ;;y2 =  1.1 * y/y(0)
; ; p2 = plot( bin_scale, y2, thick = 3, linestyle = 2, color = 'cyan', overplot = p2)

  ;;for good measure
  ;;result = ladfit(bin_scale, y/y(0))
  ;;print, 'poisson slope', result(1)
  

  ;;----------------------------------
  ;;add legend
  ;;t1 = text(1.5, 0.09, '0.02', color = 'red',/data)
  ;;t1 = text(1.5, 0.07,'0.10',color = 'coral',/data)
  ;;t1 = text(1.5, 0.055,'0.40',color = 'dark_orange',/data)
  ;;t1 = text(1.5, 0.04,'2.0',color = 'yellow',/data)
  ;;t1 = text(1.5, 0.03,'6.0',color = 'chartreuse',/data)
  ;;t1 = text(1.5, 0.025,'12',color = 'dark_green',/data)
  ;;t1 = text(1.5, 0.021,'30',color = 'cyan',/data)
  ;;t1 = text(1.5, 0.019,'100',color =  'navy',/data)


  
  TOC
end


                 ;;if aorname(a) eq 34870784 and (n lt 15) then begin
                 ;;   p1 = plot(bin_timearr, bin_corrflux, '1s', title = 'binscale'+ string(n))
                 ;;endif
                 
;              if a eq 1 and n lt 65 and n mod 3 eq 0 then begin
;                 plothist, bin_corrflux, xhist, yhist, bin = 1E-4, /noprint, /noplot
;                 ph = barplot(xhist, yhist, fill_color = 'sky blue')
;              endif
