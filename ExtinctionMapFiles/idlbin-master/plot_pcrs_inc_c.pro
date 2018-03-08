pro plot_pcrs_inc_c
   colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'blue' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
 
  ;restore, '/Users/jkrick/irac_warm/pcrsbias/pcrsbias_corr.sav'
  ;aorname = ['r44407808','r44408064','r44408320','r44408576','r44408832','r44409088','r44409344','r44409600','r44409856','r44410112']
   restore, '/Users/jkrick/irac_warm/pmap/pmap_corr.sav'
   aorname = ['r44464128','r44463872','r44463616','r44463360','r44463104','r44462848','r44462592','r44462336','r44462080','r44461824']

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

;try binning the fluxes by subarray image
nfits = 242
 bin_snaps = replicate({binsnapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname))
 

for sa = 0, n_elements(aorname) - 1 do begin
   for si = 0, nfits - 1 do begin
      idata = pmap[sa].flux[si*63:si*63 + 62]
      idataerr = pmap[sa].fluxerr[si*63:si*63 + 62]
      binned_flux = mean(idata,/nan)
      binned_fluxerr =   sqrt(total(idataerr^2))/ n_elements(idataerr)
      binned_corrflux = mean(pmap[sa].corrflux[si*63:si*63 + 62],/nan)
      binned_corrfluxerr = sqrt(total((pmap[sa].corrfluxerr[si*63:si*63 + 62])^2))/ n_elements(pmap[sa].corrfluxerr[si*63:si*63 + 62])
      binned_timearr = mean(pmap[sa].timearr[si*63:si*63 + 62])
      binned_xcen = mean(pmap[sa].xcen[si*63:si*63 + 62])
      binned_ycen = mean(pmap[sa].ycen[si*63:si*63 + 62])

      bin_snaps[sa].flux[si] = binned_flux
      bin_snaps[sa].fluxerr[si] = binned_fluxerr
      bin_snaps[sa].corrflux[si] = binned_corrflux
      bin_snaps[sa].corrfluxerr[si] = binned_corrfluxerr
      bin_snaps[sa].timearr[si] = binned_timearr
      bin_snaps[sa].xcen[si] = binned_xcen
      bin_snaps[sa].ycen[si] = binned_ycen

   endfor
   bin_snaps[sa].sclktime_0 = pmap[sa].sclktime_0

endfor

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

 z = pp_multiplot(multi_layout=[n_elements(aorname), 1], global_xtitle='X (pixels)',global_ytitle='Y (pixels)')
  for a = 0,  n_elements(aorname) - 1 do begin
     xy=z.plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.',color = 'black', thick = 2, xminor = 0, xtickinterval = 0.5)
     xy.yrange = [14.5, 15.5]
     xy.xrange = [14.5, 15.5]
  endfor

  o = pp_multiplot(multi_layout=[n_elements(aorname),1],global_xtitle='Time (hours)',global_ytitle='Flux (Jy)')
  for a = 0, n_elements(aorname) - 1 do begin
     st = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].flux,bin_snaps[a].fluxerr,  color = colorarr[a], xminor = 0, xtickinterval = 0.5, multi_index = a) ;
     st.yrange = [0.40,0.44];[1.015,1.065]
      st.xrange = [0,0.6]
    sto = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].corrflux - .01, bin_snaps[a].corrfluxerr, multi_index = a)
     sto.yrange = [0.40,0.44];[1.015,1.065]
     sto.xrange = [0,0.6]
    st1 = o.plot(bin_snaps[a].timearr, (bin_snaps[a].flux / bin_snaps[a].corrflux) -0.60, color = colorarr[a], multi_index = a)
     st1.yrange = [0.40,0.44];[1.015,1.065]
     st1.xrange = [0,0.6]
  st.xminor = 0
  endfor

;  for a = 0, n_elements(aorname) - 1 do begin;
;     if a eq 0 then begin
;        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.', thick = 2, xrange = [14.5, 15.5], yrange = [14.5, 15.5], $
;        xtitle = 'X (pixel)', ytitle = 'Y (pixel)', color = colorarr[a], aspect_ratio = 1)
;     endif

;     if a gt 0 then begin
;        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.', thick = 2, color = colorarr[a],/overplot)
;     endif
;  endfor


end

