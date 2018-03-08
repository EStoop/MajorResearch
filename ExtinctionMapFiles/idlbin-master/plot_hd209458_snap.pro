pro plot_hd209458_snap

;colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green', 'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender', 'peach_puff', 'pale_goldenrod','red']

colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'red', 'orange_red', 'dark_orange', 'gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green', 'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'red','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender', 'peach_puff', 'pale_goldenrod','green']


restore,'/Users/jkrick/irac_warm/hd209458/hd209458_snap.sav'
 aorname = ['0045188864','0045189120','0045189376','0045189632','0045189888','0045190144','0045190400','0045190656','0045190912','0045191168','0045191424','0045191680','0045191936','0045192192','0045192704','0045195264','0045192960','0045193216','0045193472','0045193984','0045193728','0045195520','0045194240','0045194496','0045194752','0045195008','0045196288','0045195776','0045197312','0045196032','0045196544','0045196800','0045197056','0045197568','0045197824','0045198080','0045192448']


 fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120827.fits', pmapdata, pmapheader
; fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/xgrid_ch2_0p1s_x4_500x500_0043_120124.fits', xgriddata, xgridheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/ygrid_ch2_0p1s_x4_500x500_0043_120124.fits', ygriddata, ygridheader

;fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/occu_ch2_0p1s_x4_500x500_0043_120124.fits',occudata, occuheader
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;try binning the fluxes by subarray image
nfits = 60
 bin_snaps = replicate({binsnapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname) )
 
sc = 0
for sa = 0, n_elements(aorname) - 1 do begin
   print, 'sa, sc, mean', sa, sc, mean(snapshots[sa].xcen), mean(snapshots[sa].ycen)
;  if sa eq 8 then begin
      ;do nothing 
      ;getting rid of this AOR
;   endif else begin
      for si = 0, nfits - 1 do begin
         idata = snapshots[sa].flux[si*63:si*63 + 62]
         idataerr = snapshots[sa].fluxerr[si*63:si*63 + 62]
         cdata = snapshots[sa].corrflux[si*63:si*63 + 62]
         cdataerr = snapshots[sa].corrfluxerr[si*63:si*63 + 62]
         binned_flux = mean(idata,/nan)
         binned_fluxerr =   sqrt(total(idataerr^2))/ n_elements(idataerr)
         binned_corrflux = mean(cdata,/nan)
         binned_corrfluxerr =  sqrt(total(cdataerr^2))/ n_elements(cdataerr)
         binned_timearr = mean(snapshots[sa].timearr[si*63:si*63 + 62])
         binned_xcen = mean(snapshots[sa].xcen[si*63:si*63 + 62])
         binned_ycen = mean(snapshots[sa].ycen[si*63:si*63 + 62])
         
         bin_snaps[sc].flux[si] = binned_flux
         bin_snaps[sc].fluxerr[si] = binned_fluxerr
         bin_snaps[sc].corrflux[si] = binned_corrflux
         bin_snaps[sc].corrfluxerr[si] = binned_corrfluxerr
         bin_snaps[sc].timearr[si] = binned_timearr
         bin_snaps[sc].xcen[si] = binned_xcen
         bin_snaps[sc].ycen[si] = binned_ycen
         
      endfor
      bin_snaps[sc].sclktime_0 = snapshots[sc].sclktime_0
      sc = sc + 1
;     endelse 
endfor

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------



;
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

  print, 'raw mean, stddev', mean(bin_snaps.flux,/nan), robust_sigma(bin_snaps.flux)
  print, 'corr mean, stddev', mean(bin_snaps.corrflux,/nan), robust_sigma(bin_snaps.corrflux)
  
;   o = pp_multiplot(multi_layout=[n_elements(aorname),1],global_xtitle='Time (hours)',global_ytitle='Flux (mJy)')
;  for a = 0,  n_elements(aorname) - 1 do begin
;     bin_snaps[a].timearr =  (bin_snaps[a].timearr - bin_snaps[a].timearr(0)) / 60./60.

;     ;print, bin_snaps[a].timearr
;     st = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].flux,bin_snaps[a].fluxerr,  color = colorarr[a], xminor = 0, xtickinterval = 0.5, multi_index = a, propagate = 0) ;
;     st.yrange = [0.485,0.51]
;     ;st.xrange = [0,0.5]
;     st.xmajor = 1
;     st.xminor = 0
;     sto = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].corrflux-0.005 , bin_snaps[a].corrfluxerr, multi_index = a, propagate = 0)
;     sto.yrange = [0.485,0.51]
;     ;sto.xrange = [0,0.5]
;     sto.xmajor = 1
;     sto.xminor = 0

;     st1 = o.plot(findgen(10) / 10., fltarr(10) + mean(bin_snaps.flux,/nan), color = 'red', multi_index = a, propagate = 0)
;     st1.yrange =[0.485,0.51]
 ;    ;st1.xrange = [0,0.5]
;     st1.xmajor = 1
;     st1.xminor = 0
;     st2 = o.plot(findgen(10) / 10., fltarr(10) + mean(bin_snaps.corrflux,/nan) - 0.005, color = 'red', multi_index = a, propagate = 0)
;     st2.yrange = [0.485,0.51]
;     ;st2.xrange = [0,0.5]
;     st2.xmajor = 1
;     st2.xminor = 0
     
;  endfor

;-------

 c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500])
; c = contour(occudata, /fill, n_levels = 7, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', title = 'occupation', aspect_ratio = 1, xrange = [0,500])

bin_snaps.xcen = 500.* (bin_snaps.xcen - 14.5)
bin_snaps.ycen = 500.* (bin_snaps.ycen - 14.5)

  for a = 0, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
;        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, xrange = [14.5, 15.5], yrange = [14.5, 15.5], xtitle = 'X (pixel)', ytitle = 'Y (pixel)', color = colorarr[a], aspect_ratio = 1)
       an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
     endif
;
     if a gt 0 then begin
        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
     endif
  endfor

  xsweet = 15.120
  ysweet = 15.085  
  box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
  box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
;  line4 = polyline(box_x, box_y, thick = 2, color = !color.black,/data)
 


  ;save, /all, filename='/Users/jkrick/irac_warm/snapshots/snaps_fit.sav'

end



  ;do some normalization
;  for a = 0, n_elements(aorname) - 1 do begin
;     snapshots[a].fluxerr = snapshots[a].fluxerr / snapshots[a].flux(1)
;     snapshots[a].flux = snapshots[a].flux / snapshots[a].flux(1)
;     snapshots[a].corrfluxerr = snapshots[a].corrfluxerr / snapshots[a].corrflux(1)
;     snapshots[a].corrflux = snapshots[a].corrflux / snapshots[a].corrflux(1)

;  endfor
