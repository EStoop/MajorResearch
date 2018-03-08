pro plot_koi69,binning = binning
  colorarr = [  'red'   , 'BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'blue' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
 
;koi69  peakup
  restore, '/Users/jkrick/irac_warm/koi69/koi69_corr.sav'
  aorname = [  'r44448512',  'r44448768']
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

if keyword_set(binning) then begin

;try binning the fluxes by subarray image
   nfits = 3549; 1521; 507.; 169.                 ; koi69 has 169 * 63 frames, must keep this constant
   nframes = 3; 7; 21.; 63.
   bin_snaps = replicate({binsnapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname))
 

   for sa = 0,  n_elements(aorname) - 1 do begin
      print, 'mean', sa, mean(snapshots[sa].xcen), mean(snapshots[sa].ycen,/nan)
      for si = 0, nfits - 1 do begin
         idata = snapshots[sa].flux[si*nframes:si*nframes + (nframes - 1)]
         idataerr = snapshots[sa].fluxerr[si*nframes:si*nframes + (nframes - 1)]
         binned_flux = mean(idata,/nan)
         binned_fluxerr =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
         binned_corrflux = mean(snapshots[sa].corrflux[si*nframes:si*nframes + (nframes - 1)],/nan)
         binned_corrfluxerr = sqrt(total((snapshots[sa].corrfluxerr[si*nframes:si*nframes + (nframes - 1)])^2))/ (n_elements(snapshots[sa].corrfluxerr[si*nframes:si*nframes + (nframes - 1)]))
         binned_timearr = mean(snapshots[sa].timearr[si*nframes:si*nframes + (nframes - 1)])
         binned_xcen = mean(snapshots[sa].xcen[si*nframes:si*nframes + (nframes - 1)])
         binned_ycen = mean(snapshots[sa].ycen[si*nframes:si*nframes + (nframes - 1)])
         
         bin_snaps[sa].flux[si] = binned_flux
         bin_snaps[sa].fluxerr[si] = binned_fluxerr
         bin_snaps[sa].corrflux[si] = binned_corrflux
         bin_snaps[sa].corrfluxerr[si] = binned_corrfluxerr
         bin_snaps[sa].timearr[si] = binned_timearr
         bin_snaps[sa].xcen[si] = binned_xcen
         bin_snaps[sa].ycen[si] = binned_ycen
         
      endfor
      bin_snaps[sa].sclktime_0 = snapshots[sa].sclktime_0
      
      
   endfor
   
endif   ;binning

if not keyword_set(binning) then begin
   bin_snaps = snapshots
endif

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

 
for a = 0, 0 do begin; n_elements(aorname) - 1 do begin
 
;getting rid of outliers
   good_pmap = where(finite(bin_snaps[a].corrflux) gt 0 and finite(bin_snaps[a].corrfluxerr) gt 0,ngood_pmap, complement=bad_pmap)
   ;print, 'bad pmap', bin_snaps[a].corrfluxerr(bad_pmap)
   xarr = bin_snaps[a].xcen[good_pmap]
   yarr = bin_snaps[a].ycen[good_pmap]
   timearr = bin_snaps[a].timearr[good_pmap]
   fluxerr = bin_snaps[a].fluxerr[good_pmap]
   flux = bin_snaps[a].flux[good_pmap]
   corrflux = bin_snaps[a].corrflux[good_pmap]
   corrfluxerr = bin_snaps[a].corrfluxerr[good_pmap]

  ;now make the plots
   
   if a eq 0 then begin
      st = plot(timearr, yarr,'1-*', yrange = [14.5,15.5],xtitle = 'Time(hours)',ytitle = 'Y pix',   color = colorarr[a])         ;
      st2 = plot(timearr, xarr,'1-*',yrange = [14.5,15.5], xtitle = 'Time(hours)',ytitle = 'X pix',  color = colorarr[a])        ;
      ;st3 = errorplot(timearr, flux,fluxerr,'1-*', xtitle = 'Time(hours)',ytitle = 'Flux',  color = colorarr[a], yrange = [0.07, 0.08] )                  ;
      ;st4 = errorplot(timearr, corrflux,corrfluxerr, '1-*', xtitle = 'Time(hours)',ytitle = 'Corrected Flux',  color = colorarr[a], yrange = [0.07, 0.08]) ; 
      st3 = plot(timearr, flux,'1-*', xtitle = 'Time(hours)',ytitle = 'Flux (Jy)',  color = colorarr[a], yrange = [0.074, 0.078] )                  ;
      st4 = plot(timearr, corrflux, '1-*', xtitle = 'Time(hours)',ytitle = 'Corrected Flux (Jy)',  color = colorarr[a], yrange = [0.074, 0.078]) ; 
      st5 = plot(xarr,  yarr, '1-*', color = colorarr[a], xrange = [14.5, 15.5], yrange = [14.5, 15.5])
      st6 = plot(timearr, flux / corrflux, '1-*', color = colorarr[a], xtitle = 'Time(hours)',ytitle = 'Correction Values')
   endif
   if a gt 0 then begin
      st.select
      st = plot(timearr, yarr,'1-*', /overplot, color = colorarr[a], /data, /current) ;
      
      st2.select
      st2 = plot(timearr, xarr,'1-*',/overplot, color = colorarr[a], /data, /current) ;
      
      st3.select
      ;st3 = errorplot(timearr, flux,fluxerr,'1-*',/overplot, color = colorarr[a], /data, /current )    
      st3 = plot(timearr, flux,'1-*',/overplot, color = colorarr[a], /data, /current )    

      st4.select
      ;st4 = errorplot(timearr, corrflux,corrfluxerr, '1-*',/overplot, color = colorarr[a], /data, /current) ; 
      st4 = plot(timearr, corrflux, '1-*',/overplot, color = colorarr[a], /data, /current) ; 

      st5.select
      st5 = plot(xarr,   yarr, '1-*', /overplot, color = colorarr[a], /data, /current)

      st6.select
      st6 = plot(timearr, flux / corrflux, '1-*', /overplot, color = colorarr[a], /data, /current)
            
   endif
  
endfor                          ;for each AOR





end


