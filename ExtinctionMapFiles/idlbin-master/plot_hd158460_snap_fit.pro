pro plot_hd158460_snap_fit
   colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' ,  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'blue' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'brown'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
 ;colorarr = [  'SILVER'   , 'SILVER'  ,'SILVER'  , 'SILVER'  ,  'SILVER'  , 'SILVER'   , 'SILVER'   , 'SILVER' , 'TEAL' ,  'teal' ,   'teal'  ,  'SILVER' , 'SILVER'   ,  'SILVER'  , 'TEAL', 'TEAL'  ,'TEAL','   PERU',   ' POWDER_BLUE' ]
 
 ; restore, '/Users/jkrick/irac_warm/snapshots/snapshots_corr.sav'
 ; aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]
  ;January 2012 snapshot take two:
;aorname = ['0045184256','0045184512','0045184768','0045185024','0045185280','0045185536','0045185792','0045186048','0045186304','0045186560','0045186816']
aorname = ['0045184256','0045184512','0045184768','0045185024','0045185280','0045185536','0045185792','0045186048','0045186304','0045186560','0045186816' ,'0045187072','0045188352','0045187584','0045187840','0045187328','0045188096','0045188608']

restore, '/Users/jkrick/irac_warm/pmap/hd158460_pmap1000.sav'


;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;try binning the fluxes by subarray image
nfits = 210;75 
 bin_snaps = replicate({binsnapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname) - 1)
 
sc = 0
for sa = 0, n_elements(aorname) - 1 do begin
   print, 'sa, sc, mean', sa, sc, mean(snapshots[sa].xcen)
   if sa eq 8 then begin
      ;do nothing 
      ;getting rid of this AOR
   endif else begin
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
     endelse 
endfor

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------



;try running this multiple times to find the statistics of success in fitting the input parameters
  nruns = 100.
  amparr = fltarr(nruns)
  amperrarr = fltarr(nruns)
  periodarr= fltarr(nruns)
  perioderrarr= fltarr(nruns)
  phasearr= fltarr(nruns)
  phaseerrarr= fltarr(nruns)
  for mc = 0., nruns - 1 do begin
 
     nsnaps = 35 ; n_elements(aorname)

;use the corrected fluxes
     combinedfluxarr = fltarr(nsnaps* n_elements(bin_snaps[0].flux))
     combinedfluxerrarr = fltarr(nsnaps * n_elements(bin_snaps[0].flux))
     xarr = fltarr(nsnaps * n_elements(bin_snaps[0].flux))
     yarr = fltarr(nsnaps * n_elements(bin_snaps[0].flux))

     ;timearr = dblarr(n_elements(aorname) * n_elements(bin_snaps[0].flux))
  ;now try 35 of them throughout a 3.5 day period
     timearr = dblarr(nsnaps * n_elements(bin_snaps[0].flux))

    ;need 35 random locations to put the
    ;snapshots between 0 and 84 hours
     random_84 = randomu(seed, nsnaps)
     random_84 = random_84 * 84.
     
    ; for a = 0, n_elements(aorname) - 1 do begin
     for a = 0, nsnaps - 1 do begin
        ;randomly pick one of the snapshots to fit to the phase
        random_15 = randomu(seed) ;make a random number between 0 and 1
        random_15 = fix(random_15 * 10.) ; now a random number between 0 and 15 ie one of the snapshot aors
        timearr[a*n_elements(bin_snaps[random_15].flux) : ((a + 1)*n_elements(bin_snaps[random_15].flux))- 1] = random_84(a) + bin_snaps[random_15].timearr

        combinedfluxarr[a*n_elements(bin_snaps[random_15].flux) : ((a + 1)*n_elements(bin_snaps[random_15].flux))- 1] = bin_snaps[random_15].corrflux ;
        combinedfluxerrarr[a*n_elements(bin_snaps[random_15].flux) : ((a + 1)*n_elements(bin_snaps[random_15].flux))- 1] = bin_snaps[random_15].corrfluxerr
        xarr[a*n_elements(bin_snaps[random_15].flux) : ((a + 1)*n_elements(bin_snaps[random_15].flux))- 1] = bin_snaps[random_15].xcen
        yarr[a*n_elements(bin_snaps[random_15].flux) : ((a + 1)*n_elements(bin_snaps[random_15].flux))- 1] = bin_snaps[random_15].ycen
        
;first step, put all the snapshot observations togehter as though they were one time sequence
     ;timearr[a*n_elements(snapshots[0].flux) : ((a + 1)*n_elements(snapshots[0].flux))- 1] = max(timearr) + 1.11E-4 + snapshots[a].timearr

;this time put them in real time order, then superpose a longer period function
     ;timearr[a*n_elements(snapshots[0].flux) : ((a + 1)*n_elements(snapshots[0].flux))- 1] = (snapshots[a].sclktime_0 + snapshots[a].timearr*60.*60.) 

; and now at random times


     endfor
  
;put back in hours for the real time timearr
;timearr =  (timearr - timearr(0)) / 60./60.

;normalize to 1 phase
     timearr = timearr / max(timearr)




;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;now add in the fake sin curve
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

;get rid of outliers in the pointing
     good = where(xarr gt 14.95 and xarr lt 15.35 and yarr gt 14.825 and yarr lt 15.2,ngood, complement = bad )
;  print, 'bad y', yarr[bad]
     xarr = xarr[good]
     yarr = yarr[good]
     timearr = timearr[good]
     combinedfluxerrarr = combinedfluxerrarr[good]
     combinedfluxarr = combinedfluxarr[good]
     
;get rid of places where there is no pmap corrections
     good_pmap = where(finite(combinedfluxarr) gt 0 and finite(combinedfluxerrarr) gt 0,ngood_pmap, complement=bad_pmap)
;     print, 'bad pmap', combinedfluxarr(bad_pmap)
     xarr = xarr[good_pmap]
     yarr = yarr[good_pmap]
     timearr = timearr[good_pmap]
     combinedfluxerrarr = combinedfluxerrarr[good_pmap]
     combinedfluxarr = combinedfluxarr[good_pmap]

;put in easy to use variables
     x = timearr
     y = combinedfluxarr
     yerr = combinedfluxerrarr
     
;the sin curve
     amplitude =  0.003          ; 0.0006
     period = 1/(2*!Pi)         ; 1 phase
     phase = 0.
     y_sin = amplitude*sin(x/period + phase ) 
     noise = randomn(seed, n_elements(x), /Double)
     yerr_sin= (noise- mean(noise))/stddev(noise)*0.0002 + 0.0002
     
     y_sum = y + y_sin
     yerr_sum = sqrt( yerr^2 + yerr_sin^2)
     
;the plots!
;     p1 = errorplot(x, y, yerr, '+2',  xtitle = 'Phase', ytitle = 'Flux (Jy) ',  xrange = [0,1], yrange = [1.02, 1.08]) ;, errorbar_capsize = .001)
;     p3 = errorplot(x, y_sum,yerr_sum, '+b2',/overplot)


;can I recover the sin curve? Maybe?
;don't need to fit for all the other stuff

     print, 'n_elements', n_elements(x), n_elements(y_sum), n_elements(yerr_sum)
     print, 'mean', mean(x), mean(y_sum), mean(yerr), mean(yerr_sin)

     start = [amplitude, mean(y_sum), period,phase]
  ;limit the range of the sin cuve phase
     parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start))
     parinfo[3].fixed = 1 
     parinfo[2].fixed = 1
     parinfo[1].limited[0] = 1
     parinfo[1].limits[0] = 0.0
     pa= MPFITFUN('sin_func',x,y_sum, yerr_sum, start,PARINFO =parinfo , PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg,/quiet) 
     fitx = findgen(100) / 100.
;     p4 = plot( fitx, pa(0)*sin(fitx/pa(2) + pa(3)) + pa(1),  '6r1--',/overplot) ;fitted curve
;     p4 = plot( fitx, amplitude*sin(fitx/period +phase)+ mean(y_sum),  'r1',/overplot)   ;input curve
  
     print, 'result', pa
     print, 'status', status
     print, errmsg
     achi = achi / adof
     print, 'reduced chi squared', achi
     print, 'input amplitude, period, phase', amplitude, period*2.*!pi, phase/(2.*!pi)
     print, 'fitted amplitude', pa[0], spa[0]
     ;print, 'fitted period (hours)', pa[2]*2.*!pi, spa[2]*2.*!pi
     ;print, 'fitted phase [0,1]', pa[3]/(2.*!pi), spa[3]/(2.*!pi)
     
     ;keep track of the results
     amparr[mc] = pa[0]
     amperrarr[mc] = spa[0]
     periodarr[mc] = pa[2]
     perioderrarr[mc] = spa[2]
     phasearr[mc] = pa[3]
     phaseerrarr[mc] = spa[3]

 ;plot residuals
     sub = y_sum - ( pa(0)*sin(x/pa(2) + pa(3)) + pa(1))
;     f2 = plot(timearr, sub+1.03, '6g1.',/overplot)


;fit residuals to know what the errors look like.
  ;;   plothist, sub, xhist, yhist,bin=0.0010, /noprint,/noplot
;  r0 = barplot(xhist, yhist, fill_color = 'green', xrange = [-0.05, 0.05], xtitle = 'Residuals', ytitle = 'Number')

  ;;   noise = fltarr(n_elements(yhist)) + 1.0 ;  
  ;;   by = where(yhist lt 100.)
  ;;   noise(by) = 100.
  ;;   start = [0.0001, 0.01, 10000.]
  ;;   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) 
;  r1 = plot( xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.),  '3',/overplot)
  
;  t_sig = text(0.02, 2000, string(result(1)),'3',/data)

  endfor  ;for the mc iterations

;what are the statistics
  success_amp = where(amparr gt (amplitude - 0.2*amplitude) and amparr lt (amplitude + 0.2*amplitude), nsuccessamp) ;within 20% of the real answer
  ;print, 'nsuccessamp', nsuccessamp
  print, (nsuccessamp/ mc) * 100., '% of the runs measured amplitude within 20% of the input value'

;  success_per = where(periodarr gt (period - 0.2*period) and periodarr lt (period + 0.2*period), nsuccessper) ;within 20% of the real answer
  ;print, 'nsuccessaper', nsuccessper
;  print, (nsuccessper/ mc) * 100., '% of the runs measured period within 20% of the input value'

;  success_phase = where(phasearr gt (phase - 0.2*phase) and phasearr lt (phase + 0.2*phase), nsuccessphase) ;within 20% of the real answer
  ;print, 'nsuccessphase', nsuccessphase
;  print, (nsuccessphase/ mc) * 100., '% of the runs measured phase within 20% of the input value'


;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

;  z = pp_multiplot(multi_layout=[n_elements(aorname), 1], global_xtitle='X (pixels)',global_ytitle='Y (pixels)')
;  for a = 0,  n_elements(aorname) - 1 do begin
;     xy=z.plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1.',color = 'black', thick = 2, xminor = 0, xtickinterval = 0.5)
;     xy.yrange = [14.5, 15.5]
;     xy.xrange = [14.5, 15.5]
;  endfor
;

  print, 'raw mean, stddev', mean(bin_snaps.flux,/nan), robust_sigma(bin_snaps.flux)
  print, 'corr mean, stddev', mean(bin_snaps.corrflux,/nan), robust_sigma(bin_snaps.corrflux)
  
  o = pp_multiplot(multi_layout=[n_elements(aorname),1],global_xtitle='Time (hours)',global_ytitle='Flux (mJy)')
  for a = 0,  n_elements(aorname) - 2 do begin
     st = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].flux,bin_snaps[a].fluxerr,  color = colorarr[a], xminor = 0, xtickinterval = 0.5, multi_index = a) ;
     st.yrange = [1.03,1.06]    ;[1.015,1.065]
     st.xrange = [0,0.5]
     st.xmajor = 1
     st.xminor = 0
     sto = o.errorplot(bin_snaps[a].timearr, bin_snaps[a].corrflux-0.005 , bin_snaps[a].corrfluxerr, multi_index = a)
     sto.yrange = [1.03,1.06]   ;[1.015,1.065]
     sto.xrange = [0,0.5]
     sto.xmajor = 1
     sto.xminor = 0
                                ; st1 = o.plot(bin_snaps[a].timearr, (bin_snaps[a].flux / bin_snaps[a].corrflux) + 0.02, color = colorarr[a], multi_index = a)
   ;  st1.yrange = [1.02,1.06];[1.015,1.065]
   ;  st1.xrange = [0,0.6]


     st1 = o.plot(findgen(10) / 10., fltarr(10) + mean(bin_snaps.flux,/nan), color = 'red', multi_index = a)
     st1.yrange = [1.03,1.06]   ;[1.015,1.065]
     st1.xrange = [0,0.5]
     st1.xmajor = 1
     st1.xminor = 0
     st2 = o.plot(findgen(10) / 10., fltarr(10) + mean(bin_snaps.corrflux,/nan) - 0.005, color = 'red', multi_index = a)
     st2.yrange = [1.03,1.06]   ;[1.015,1.065]
     st2.xrange = [0,0.5]
     st2.xmajor = 1
     st2.xminor = 0
     
  endfor

  for a = 0, n_elements(aorname) - 2 do begin
     if a eq 0 then begin
        an = plot(bin_snaps[a].xcen, bin_snaps[a].ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, xrange = [14.5, 15.5], yrange = [14.5, 15.5], xtitle = 'X (pixel)', ytitle = 'Y (pixel)', color = colorarr[a], aspect_ratio = 1)
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
  line4 = polyline(box_x, box_y, thick = 2, color = !color.black,/data)
 


  save, /all, filename='/Users/jkrick/irac_warm/snapshots/snaps_fit_1000pmap.sav'

end



  ;do some normalization
;  for a = 0, n_elements(aorname) - 1 do begin
;     snapshots[a].fluxerr = snapshots[a].fluxerr / snapshots[a].flux(1)
;     snapshots[a].flux = snapshots[a].flux / snapshots[a].flux(1)
;     snapshots[a].corrfluxerr = snapshots[a].corrfluxerr / snapshots[a].corrflux(1)
;     snapshots[a].corrflux = snapshots[a].corrflux / snapshots[a].corrflux(1)

;  endfor
