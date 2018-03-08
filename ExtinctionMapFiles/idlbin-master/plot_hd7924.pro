pro plot_hd7924,binning = binning
  colorarr = ['black',   'red'   , 'BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'blue' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
 
;koi69  peakup
  restore, '/Users/jkrick/irac_warm/snapshots/hd7924/hd7924_corr.sav'
  aorname = ['r44605184']
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
  
  if keyword_set(binning) then begin
     
;try binning the fluxes by subarray image
     nfits = 17.                ; 85;170.;340.;680.;850;1700.;3400                
     nframes =12600.            ; 2520.;1260.;630.;315.;252.;126.;63.
     bin_snaps = replicate({binsnapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ',back:fltarr(nfits)},n_elements(aorname))
     
     
     for si = 0., nfits - 1 do begin
                                ;print, 'working on si', si
        idata = snapshots[0].flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = snapshots[0].fluxerr[si*nframes:si*nframes + (nframes - 1)]
        iback = snapshots[0].back[si*nframes:si*nframes + (nframes - 1)]
        binned_flux = mean(idata,/nan)
        binned_fluxerr =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        binned_back = mean(iback,/nan)
        binned_corrflux = mean(snapshots[0].corrflux[si*nframes:si*nframes + (nframes - 1)],/nan)
        binned_corrfluxerr = sqrt(total((snapshots[0].corrfluxerr[si*nframes:si*nframes + (nframes - 1)])^2))/ (n_elements(snapshots[0].corrfluxerr[si*nframes:si*nframes + (nframes - 1)]))
        binned_timearr = mean(snapshots[0].timearr[si*nframes:si*nframes + (nframes - 1)])
        binned_xcen = mean(snapshots[0].xcen[si*nframes:si*nframes + (nframes - 1)])
        binned_ycen = mean(snapshots[0].ycen[si*nframes:si*nframes + (nframes - 1)])
        
        bin_snaps[0].flux[si] = binned_flux
        bin_snaps[0].fluxerr[si] = binned_fluxerr
        bin_snaps[0].corrflux[si] = binned_corrflux
        bin_snaps[0].corrfluxerr[si] = binned_corrfluxerr
        bin_snaps[0].timearr[si] = binned_timearr
        bin_snaps[0].xcen[si] = binned_xcen
        bin_snaps[0].ycen[si] = binned_ycen
        bin_snaps[0].back[si] = binned_back
        
     endfor                     ;for each fits image
     bin_snaps[0].sclktime_0 = snapshots[0].sclktime_0
 
  endif                         ;binning
  
  if not keyword_set(binning) then begin
     bin_snaps = snapshots
  endif
  
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

 
  a = 0                                                                                   ;only 1 AOR 
;getting rid of outliers
  good_pmap = where(finite(bin_snaps[a].corrflux) gt 0 ,ngood_pmap, complement=bad_pmap) ;and finite(bin_snaps[a].corrfluxerr) gt 0
   ;print, 'bad pmap', bin_snaps[a].corrfluxerr(bad_pmap)
  xarr = bin_snaps[a].xcen[good_pmap]
  yarr = bin_snaps[a].ycen[good_pmap]
  timearr = bin_snaps[a].timearr[good_pmap]
  fluxerr = bin_snaps[a].fluxerr[good_pmap]
  flux = bin_snaps[a].flux[good_pmap]
  corrflux = bin_snaps[a].corrflux[good_pmap]
  corrfluxerr = bin_snaps[a].corrfluxerr[good_pmap]
  backarr = bin_snaps[a].back[good_pmap]
  
   ;now make the plots
   
      ;-----------------
      ;y vs. time
;      st = plot(timearr, yarr,'1-o', yrange = [14.9,15.1],xtitle = 'Time(hours)',ytitle = 'Y pix',   color = colorarr[a], sym_filled = 1, sym_size = 0.4, xrange = [2,4])         ;
;      pl = polyline([2.7,2.7], [0,100], '--', color= 'black',/data)
;      pl = polyline([3.1,3.1], [0,100], '--', color= 'black',/data)

      ;-----------------
      ;x vs time
;      st2 = plot(timearr, xarr,'1-o',yrange = [15.0,15.2], xtitle = 'Time(hours)',ytitle = 'X pix',  color = colorarr[a], sym_filled = 1, sym_size = 0.4, xrange = [2,4])        ;
;      pl = polyline([2.7,2.7], [0,100], '--', color= 'black',/data)
;      pl = polyline([3.1,3.1], [0,100], '--', color= 'black',/data)

      ;-----------------
      ;flux vs time
;      st3 = plot(timearr, flux,'1-', xtitle = 'Time(hours)',ytitle = 'Flux (Jy)',  title = 'Pmap', color = 'light gray' , yrange = [1.435, 1.465]); , xrange = [2, 4]) ;
;      st3 = plot(timearr, flux,'1o',  color = 'black' , sym_filled = 1, sym_size = 0.3 , /overplot)        ; , xrange = [2, 4])                  ;
;      st4 = plot(timearr, corrflux +0.01, '1-',  color = 'light gray',/overplot) ; 
;      st4 = plot(timearr, corrflux +0.01, '1o',  color = 'red', sym_filled = 1, sym_size = 0.3 ,/overplot) ; 
;      pl = polyline([2.7,2.7], [0,100], '--', color= 'black',/data)
;      pl = polyline([3.1,3.1], [0,100], '--', color= 'black',/data)


      ;-----------------
                                ;Normalized flux vs time
                                ; st6 = plot(timearr, corrflux/corrflux(0), '1-o',  color = 'red', sym_filled = 1, sym_size = 0.4,NAME = 'HD7924', yrange = [0.998, 1.0023], xrange = [0,8],xtitle = 'Time (hours)', ytitle = 'Normalized Flux ', title = 'PMAP')

  
  st6 = errorplot(timearr, corrflux/median(corrflux), corrfluxerr/median(corrflux), '1-o', color = 'red', sym_filled = 1, sym_size = 0.4,NAME = 'pmap', yrange = [0.998, 1.002], xrange = [0,8],xtitle = 'Time (hours)', ytitle = 'Normalized Flux ', title = 'HD7924')
;this is the uncorrected flux
  p1 = errorplot(timearr, flux/median(flux), fluxerr/median(flux), '1-o', color = 'black', sym_size = 0.4, sym_filled = 1.,/overplot, NAME = 'uncorrected') ;, xrange = [2,4]

     ;-----------------
      ;normalized flux  vs time including background
;     st2 = plot(timearr, backarr/mean(backarr,/nan), '1-', color = 'light gray',  ytitle = 'Normalized  flux', xtitle = 'Time(Hours)',yrange = [0.9, 1.1])
;     st5 = plot(timearr, backarr/mean(backarr,/nan), '1o', color = 'blue', sym_filled = 1, sym_size = 0.4 ,NAME = 'Background',/overplot)
;      st6 = plot(timearr, corrflux/corrflux(0), '1o',  color = 'red',/overplot, sym_filled = 1, sym_size = 0.4,NAME = 'HD7924')
;      l = legend(Target = [st5,st6], position = [5, 1.09], /data)

     ;-----------------
      ;correction values vs. time
    ;;  st6 = plot(timearr, (flux / corrflux) + 0.465, '6b1.', /overplot)
      ;st6 = plot(timearr, flux / corrflux, '1-*', color = colorarr[a], xtitle = 'Time(hours)',ytitle = 'Correction Values')
     ; st5 = plot(xarr,  yarr, '1-*', color = colorarr[a], xrange = [14.5, 15.5], yrange = [14.5, 15.5])
 


;need to print out filenanes where the blip is

;good = where(bin_snaps[0].corrflux+0.01 lt 1.44)
;print,'frames', good
;print, 'ncorr', n_elements(bin_snaps[0].corrflux)
      
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
;now work on the data fitting with a polynomial for comparison:

   restore,  '/Users/jkrick/irac_warm/snapshots/hd7924/hd7924_binnedsub.sav'
   pd = plot(timearr, (binned_sub + median(flux))/ median(binned_sub+median(flux)), '1-o', color = 'blue', sym_size = 0.4, sym_filled = 1.,/overplot, NAME = 'data fit')
   l = legend(Target = [p1,st6,pd], position = [5, 1.0016], /data)

end
