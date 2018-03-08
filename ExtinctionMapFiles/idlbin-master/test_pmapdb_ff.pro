pro test_pmapdb_ff, planetname, run_photometry = run_photometry,make_plots = make_plots
 t1 = systime(1)

;use HD20003 r48408064 for testing = aor 0
;photometry already run with xfwhm and yfwhm
;  planetname = 'HD7924b';'HD20003'
  chname = '2' & apradius = 2.25 & bin_level = 63L*2 & nn = 50
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/') 
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)
  pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140805.sav'
;'/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140716.sav'

  if chname eq '2' then occ_filename =  '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120827_occthresh.fits'$
  else occ_filename = '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch1_500x500_0043_120828_occthresh.fits'
  fits_read,occ_filename, occdata, occheader

;do a few up-front tests on the pmap database
  restore, pmapfile
  cgHistoplot, f_pmap,/nan , thick = 2, xtitle = 'Flux', ytitle = 'Number',  $
               datacolorname = 'cyan', title = 'PMAP databases',/outline, $
               HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize, xrange = [0.41, 0.44]
  binCenters = loc + (binsize / 2.0)
  yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
  fwhm = 2*SQRT(2*ALOG(2))*coeff(2)
  cgtext, 0.411, 13000, '140805 fwhm: ' + strnsignif(fwhm, 2), color = 'cyan'
  cgPlot, binCenters, yfit, COLOR='cyan', THICK=1, /OVERPLOT

  pmapfileold = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140716.sav'
  restore, pmapfileold
  cgHistoplot, f_pmap,/nan , thick = 2,datacolorname = 'blue', $
               HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize, /outline,/oplot
  binCenters = loc + (binsize / 2.0)
  yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
  fwhm = 2*SQRT(2*ALOG(2))*coeff(2)
  cgtext, 0.411, 15000, '140716 fwhm: ' + strnsignif(fwhm, 2), color = 'blue'
  cgPlot, binCenters, yfit, COLOR='blue', THICK=1, /OVERPLOT
;--------------------------------------------------------------------

  if keyword_set(run_photometry) then begin
;--------------------------------------------------------------------
;2) pmap_correct with x, y, np
     corrfluxff = pmap_correct((planethash[aorname(0),'xcen']),(planethash[aorname(0),'ycen']),planethash[aorname(0),'flux'],$
                              2,planethash[aorname(0),'xfwhm'],planethash[aorname(0),'yfwhm'],NP = planethash[aorname(0),'npcentroids'],$
                              FUNC=planethash[aorname(0),'fluxerr'],CORR_UNC=corr_unc,$
                              DATAFILE=pmapfile,NNEAREST=nn)
     save, corrfluxff, filename = '/Users/jkrick/irac_warm/pcrs_planets/' + planetname + '/corrflux_ff.sav'
     print, 'finished pmap_correct x, y, np'
     
     restore, '/Users/jkrick/irac_warm/pcrs_planets/'+planetname+'/corrflux2.sav'

;     corrflux2 = pmap_correct((planethash[aorname(0),'xcen']),(planethash[aorname(0),'ycen']),planethash[aorname(0),'flux'],$
;                              2,planethash[aorname(0),'xfwhm'],planethash[aorname(0),'yfwhm'],NP = planethash[aorname(0),'npcentroids'],$
;                              FUNC=planethash[aorname(0),'fluxerr'],CORR_UNC=corr_unc,$
;                              DATAFILE=pmapfile,NNEAREST=nn)
;     corrfluxmult = corrflux2
;     save, corrflux2, filename = '/Users/jkrick/irac_warm/pcrs_planets/' + planetname + '/corrflux2.sav'
;     print, 'finished pmap_correct x, y, np'

;--------------------------------------------------------------------
;4) nearest_neighbors on the data itself with x, y
;     pixphasecorr_noisepix, planetname, nn, apradius,chname, /use_np
;     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_np.sav',/remove_all)
;     restore, savename
;     nnflux2 = flux
     
;--------------------------------------------------------------------
  endif else begin              ;keyword_set(run_photometry)
     restore, '/Users/jkrick/irac_warm/pcrs_planets/'+planetname+'/corrflux2.sav'
     restore, '/Users/jkrick/irac_warm/pcrs_planets/'+planetname+'/corrflux_ff.sav'
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
;     restore, savename

;     nnflux1 = flux
;     nnflux2 = flux_np
;     nnflux3 = flux_fwhm
  endelse                       ;not running the photometry
  
 
;binning
  timearr =  (planethash[aorname(0),'timearr'] - (planethash[aorname(0),'timearr'])(0))
  fluxarr = planethash[aorname(0),'flux'] 
  numberarr = findgen(n_elements(fluxarr))
  h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  
  bin_time = dblarr(n_elements(h))
  bin_flux = bin_time
  bin_corrflux2= bin_time
  bin_corrfluxff= bin_time
;  bin_nn1 = bin_time
;  bin_nn2 = bin_time
;  bin_nn3 = bin_time

  c = 0
  for j = 0L, n_elements(h) - 1 do begin
     if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        bin_time[c] = median(timearr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, fluxarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_flux[c] = meanx     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux2[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrfluxff[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrfluxff[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;        meanclip_jk, nnflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;        bin_nn1[c] = meanx      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;        meanclip_jk, nnflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;        bin_nn2[c] = meanx      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;        meanclip_jk, nnflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;        bin_nn3[c] = meanx      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        c = c + 1

     endif
  endfor
  bin_time = bin_time[0:c-1]
  bin_flux = bin_flux[0:c-1]
  bin_corrflux2 = bin_corrflux2[0:c-1]
  bin_corrfluxff= bin_corrfluxff[0:c-1]
;  bin_nn1 = bin_nn1[0:c-1]
;  bin_nn2 = bin_nn2[0:c-1]
;  bin_nn3 = bin_nn3[0:c-1]
;  bin_corrfluxadd = bin_corrfluxadd[0:c-1]
;  bin_corrfluxmult = bin_corrfluxmult[0:c-1]
  print, 'finished binning'
;--------------------------------------------------------------------
;plot the raw fluxes
;make histograms to compare widths for the different methods
;    print, 'n_elements nn', n_elements(where(finite(fluxarr))), n_elements(where(finite(nnflux1))), n_elements(where(finite(nnflux2))), n_elements(where(finite(nnflux3)))
    print, 'n_elements corr', n_elements(where(finite(fluxarr))),n_elements(where(finite(corrflux2))), n_elements(where(finite(corrfluxff)))

;    cgPS_Open, '/Users/jkrick/irac_warm/pcrs_planets/'+planetname+'/test_pmap_hist.ps'

    cgHistoplot, fluxarr/ mean(fluxarr,/nan),/nan , thick = 2, xtitle = 'Flux', ytitle = 'Number',  $
                 datacolorname = 'black', title = planetname,/outline,xrange = [0.985, 1.015] , $
                 HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize ;, yrange = [1E-5, 1] , peak_height = 1.0
    binCenters = loc + (binsize / 2.0)
    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
    fwhmraw = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='black', THICK=3, /OVERPLOT
    cgtext, 0.985, 3100, string(fwhmraw), color = 'black'


;    cgHistoplot, nnflux1/ mean(nnflux1,/nan),/nan, thick = 2 , datacolorname = 'gray',/oplot,/outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
;    binCenters = loc + (binsize / 2.0)
;    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
;    fwhmnn1 = 2*SQRT(2*ALOG(2))*coeff(2)
; ;   cgPlot, binCenters, yfit, COLOR='gray', THICK=3, /OVERPLOT
;    cgtext, 0.985, 2900, string(fwhmnn1), color = 'gray'

;    cgHistoplot, nnflux2/ mean(nnflux2,/nan),/nan , thick = 2 ,datacolorname = 'red',/oplot,/outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
;    binCenters = loc + (binsize / 2.0)
;    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
;    fwhmnn2 = 2*SQRT(2*ALOG(2))*coeff(2)
;;    cgPlot, binCenters, yfit, COLOR='red', THICK=3, /OVERPLOT
;    cgtext, 0.985, 2700, string(fwhmnn2), color = 'red'

;    cgHistoplot, nnflux3/ mean(nnflux3,/nan),/nan, thick = 2 , datacolorname = 'salmon',/oplot,/outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
;    binCenters = loc + (binsize / 2.0)
;    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
;    fwhmnn3 = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='salmon', THICK=3, /OVERPLOT
;    cgtext, 0.985, 2500, string(fwhmnn3), color = 'salmon'


    cgHistoplot, corrflux2/ mean(corrflux2,/nan),/nan, thick = 2 ,datacolorname = 'blue',/oplot,/outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
    binCenters = loc + (binsize / 2.0)
    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
    fwhmc2 = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='blue', THICK=3, /OVERPLOT
    cgtext,0.985, 2100,string(fwhmc2), color = 'blue'

    cgHistoplot, corrfluxff/ mean(corrfluxff,/nan),/nan, thick = 2 , datacolorname = 'cyan',/oplot,/outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
    binCenters = loc + (binsize / 2.0)
    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
    fwhmc3 =  2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='cyan', THICK=3, /OVERPLOT
    cgtext, 0.985, 1900,string(fwhmc3), color = 'cyan'



;   cgPS_Close
  

  if keyword_set(make_plots) then begin
;pmap_correct
     p = plot(bin_time/60./60., bin_flux / mean(bin_flux), xtitle = 'Time', ytitle = 'Normalized flux', '1s', title = planetname +  ' pmap correct', sym_size = 0.2,sym_filled = 1,color = 'black', name = 'raw'+ string(fwhmraw))
     p2= plot(bin_time/60./60., (bin_corrflux2 / mean(bin_corrflux2,/nan)) + 0.005,  '1s', sym_size = 0.3,sym_filled = 1,color = 'blue', overplot = p, name = 'X Y NP'+ string(fwhmc2))
     p3= plot(bin_time/60./60., (bin_corrfluxff / mean(bin_corrfluxff,/nan)) + 0.005,  '1s', sym_size = 0.3,sym_filled = 1,color = 'cyan', overplot = p, name = 'X Y NP FF'+ string(fwhmc3))
     l = legend(target = [p,  p2, p3], position = [4, 0.999], /data, /auto_text_color)
     
 
;pixphasecorr_noisepix = nn
;     p = plot(bin_time/60./60., bin_flux / mean(bin_flux), xtitle = 'Time', ytitle = 'Normalized flux', title = planetname + ' nearest neighbor', '1s', sym_size = 0.2,sym_filled = 1,color = 'black', name = 'raw' + string(fwhmraw))
;     p4= plot(bin_time/60./60., (bin_nn1 / mean(bin_nn1,/nan)) + 0.005,  '1s', sym_size = 0.3,sym_filled = 1,color = 'gray', overplot = p, name = 'X Y'+ string(fwhmnn1))
;     p5= plot(bin_time/60./60., (bin_nn2 / mean(bin_nn2,/nan)) + 0.005,  '1s', sym_size = 0.3,sym_filled = 1,color = 'red', overplot = p, name = 'X Y NP'+ string(fwhmnn2))
;     p6= plot(bin_time/60./60., (bin_nn3 / mean(bin_nn3,/nan)) + 0.005,  '1s', sym_size = 0.3,sym_filled = 1,color = 'salmon', overplot = p, name = 'X Y XYFWHM'+ string(fwhmnn3))
;     l = legend(target = [p, p4, p5, p6], position = [6, 0.998], /data, /auto_text_color)
     
 endif   ;make_plots
  jumpend: print, 'done'
 end

