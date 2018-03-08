pro plot_hd158460_stare_v2,binning = binning
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'blue' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
 
;  restore, '/Users/jkrick/irac_warm/snapshots/stare_corr.sav'
;  aorname = ['r42506496','r42051584']
;koi69  peakup
  restore, '/Users/jkrick/irac_warm/koi69/koi69_corr.sav'
  aorname = [ 'r44448000', 'r44448256',  'r44448512',  'r44448768']
  yrange = [0.06, 0.08]
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

if keyword_set(binning) then begin

;try binning the fluxes by subarray image
   nfits = 1330.                ; staring hd158460
   nfits = 169.                 ; koi69

   bin_snaps = replicate({binsnapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname))
 

   for sa = 0,  n_elements(aorname) - 1 do begin
      print, 'mean', sa, mean(snapshots[sa].xcen), mean(snapshots[sa].ycen,/nan)
      for si = 0, nfits - 1 do begin
         idata = snapshots[sa].flux[si*63.:si*63. + 62]
         idataerr = snapshots[sa].fluxerr[si*63.:si*63. + 62]
         binned_flux = mean(idata,/nan)
         binned_fluxerr =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
         binned_corrflux = mean(snapshots[sa].corrflux[si*63.:si*63. + 62],/nan)
         binned_corrfluxerr = sqrt(total((snapshots[sa].corrfluxerr[si*63.:si*63. + 62])^2))/ (n_elements(snapshots[sa].corrfluxerr[si*63.:si*63. + 62]))
         binned_timearr = mean(snapshots[sa].timearr[si*63.:si*63. + 62])
         binned_xcen = mean(snapshots[sa].xcen[si*63.:si*63. + 62])
         binned_ycen = mean(snapshots[sa].ycen[si*63.:si*63. + 62])
         
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

  ;z = pp_multiplot(multi_layout=[n_elements(aorname), 1], global_xtitle='X (pixels)',global_ytitle='Y (pixels)')
  z = pp_multiplot(multi_layout=[2, 1], global_xtitle='X (pixels)',global_ytitle='Y (pixels)')
  for a = 2,  n_elements(aorname) - 1 do begin
     xy=z.plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.',color = 'black', thick = 2, xminor = 0, xtickinterval = 0.5)
     xy.yrange = [14.5, 15.5]
     xy.xrange = [14.5, 15.5]
  endfor

  ;o = pp_multiplot(multi_layout=[n_elements(aorname),1],global_xtitle='Time (hours)',global_ytitle='Flux (Jy)')
  o = pp_multiplot(multi_layout=[2,1],global_xtitle='Time (hours)',global_ytitle='Flux (Jy)')
  for a = 2,  n_elements(aorname) - 1 do begin
     ;st = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].flux,bin_snaps[a].fluxerr,  color = colorarr[a], xminor = 0, xtickinterval = 0.5, multi_index = a) ;
     st = o.plot(bin_snaps[a].timearr, bin_snaps[a].flux, color = colorarr[a], xminor = 0, xtickinterval = 0.5, multi_index = a) ;
     st.yrange = yrange ; [1.02,1.06];[1.015,1.065]
      st.xrange = [0,7.]
     st.xminor = 0
    ;sto = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].corrflux - .01, bin_snaps[a].corrfluxerr, multi_index = a)
    sto = o.plot(bin_snaps[a].timearr, bin_snaps[a].corrflux - .01, multi_index = a)
     sto.yrange = yrange; [1.02,1.06];[1.015,1.065]
     sto.xrange = [0,7.]
    st1 = o.plot(bin_snaps[a].timearr, (bin_snaps[a].flux / bin_snaps[a].corrflux) -0.93, color = colorarr[a], multi_index = a)
    print, 'mean corr', mean( (bin_snaps[a].flux / bin_snaps[a].corrflux),/nan)
     st1.yrange = yrange;[1.02,1.06];[1.015,1.065]
     st1.xrange = [0,7.]


  endfor

  for a = 2, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.', thick = 2, xrange = [14.5, 15.5], yrange = [14.5, 15.5], xtitle = 'X (pixel)', ytitle = 'Y (pixel)', color = colorarr[a]);, aspect_ratio = 1)
     endif

     if a gt 0 then begin
        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.', thick = 2, color = colorarr[a],/overplot)
     endif
  endfor


;for a = 0, n_elements(aorname) - 1 do begin
 
;getting rid of outliers
;   good_pmap = where(finite(bin_snaps[a].corrflux) gt 0 and finite(bin_snaps[a].corrfluxerr) gt 0,ngood_pmap, complement=bad_pmap)
;   print, 'bad pmap', bin_snaps[a].corrfluxerr(bad_pmap)
;   xarr = bin_snaps[a].xcen[good_pmap]
;   yarr = bin_snaps[a].xcen[good_pmap]
;   timearr = bin_snaps[a].timearr[good_pmap]
;   fluxerr = bin_snaps[a].fluxerr[good_pmap]
;   flux = bin_snaps[a].flux[good_pmap]
;   corrflux = bin_snaps[a].corrflux[good_pmap]
;   corrfluxerr = bin_snaps[a].corrfluxerr[good_pmap]

   ;now make the plots
;   st = plot(timearr, yarr,'6r1.', yrange = [14.5,15.5],xtitle = 'Time(hours)',ytitle = 'Y pix',  title = aorname(a)) ;
;   st2 = plot(timearr, xarr,'6r1.',yrange = [14.5,15.5], xtitle = 'Time(hours)',ytitle = 'X pix',  title = aorname(a)) ;
;   st3 = errorplot(timearr, flux,fluxerr,'6r1.', xtitle = 'Time(hours)',ytitle = 'Flux',  title =aorname(a))          ;
;   st4 = errorplot(timearr, corrflux,corrfluxerr, '6r1.', xtitle = 'Time(hours)',ytitle = 'Corrected Flux',  title = aorname(a)) ;
   
;   if a eq 0 then begin

;      st5 = plot(bin_snaps[a].xcen,   bin_snaps[a].ycen, '6.', color = colorarr[a])

;   endif
;   if a gt 0 then begin
;      st5.select
;      st5 = plot(bin_snaps[a].xcen,   bin_snaps[a].ycen, '6.', /overplot, color = colorarr[a], /data, /current)
;   endif
   
;endfor                          ;for each AOR





end


