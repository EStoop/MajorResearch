pro test_pmap_correct, planetname, startaor, stopaor, run_nn = run_nn, run_pmap = run_pmap,plot_pmap = plot_pmap, plot_nn = plot_nn, pmapff = pmapff, halfweights = halfweights
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
  print, 'aorname', aorname(startaor), aorname(stopaor)
  ;4xbinning, including all 64 subframes ; this is the standard
  pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r3p00_s3_7_0p1s_x4_140314.sav'
  pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140716.sav'
  ;4xbinning, removing first four subframes
  if keyword_set(pmapff) then pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140805.sav'

  if chname eq '2' then occ_filename =  '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120827_occthresh.fits'$
  else occ_filename = '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch1_500x500_0043_120828_occthresh.fits'
  fits_read,occ_filename, occdata, occheader

  if keyword_set(run_pmap) then begin
;--------------------------------------------------------------------
;3) pmap_correct with x, y, xfwhm, yfwhm
     corrflux3 = pmap_correct((planethash[aorname(startaor),'xcen']),(planethash[aorname(startaor),'ycen']),$
                              planethash[aorname(startaor),'flux'],2,planethash[aorname(startaor),'xfwhm'],$
                              planethash[aorname(startaor),'yfwhm'],$
                              FUNC=planethash[aorname(startaor),'fluxerr'],CORR_UNC=corr_unc,$
                              DATAFILE=pmapfile,NNEAREST=nn)
     savefile =strcompress(dirname+ '/corrflux3.sav',/remove_all)
     if keyword_set(pmapff) then savefile = strcompress(dirname + '/corrflux3_noff.sav',/remove_all)
     save, corrflux3, filename =savefile
     print, 'finished pmap_correct x, y, xyfwhm'
;--------------------------------------------------------------------
     
;1) pmap_correct with x, y
     corrflux1= pmap_correct((planethash[aorname(startaor),'xcen']),(planethash[aorname(startaor),'ycen']),$
                             planethash[aorname(startaor),'flux'],2,planethash[aorname(startaor),'xfwhm'],$
                             planethash[aorname(startaor),'yfwhm'],/POSITION_ONLY,$
                             FUNC=planethash[aorname(startaor),'fluxerr'],CORR_UNC=corr_unc,$
                             DATAFILE=pmapfile,NNEAREST=nn) 
     savefile =strcompress(dirname + '/corrflux1.sav',/remove_all)
     if keyword_set(pmapff) then savefile = strcompress(dirname + '/corrflux1_noff.sav',/remove_all)
     save, corrflux1, filename =savefile
     print, 'finished pmap_correct x, y'
;--------------------------------------------------------------------
;2) pmap_correct with x, y, np
     corrflux2 = pmap_correct((planethash[aorname(startaor),'xcen']),(planethash[aorname(startaor),'ycen']),$
                              planethash[aorname(startaor),'flux'],2,planethash[aorname(startaor),'xfwhm'],$
                              planethash[aorname(startaor),'yfwhm'],NP = planethash[aorname(startaor),'npcentroids'],$
                              FUNC=planethash[aorname(startaor),'fluxerr'],CORR_UNC=corr_unc,$
                              DATAFILE=pmapfile,NNEAREST=nn)
     savefile =strcompress(dirname + '/corrflux2.sav',/remove_all)
     if keyword_set(pmapff) then savefile = strcompress(dirname + '/corrflux2_noff.sav',/remove_all)
     save, corrflux2, filename =savefile
     print, 'finished pmap_correct x, y, np'
  endif

;--------------------------------------------------------------------
  if keyword_set(run_nn) then begin
;4) nearest_neighbors on the data itself with x, y
     pixphasecorr_noisepix, planetname, nn, apradius,chname,startaor, stopaor, /xyonly
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_xy.sav',/remove_all)
     restore, savename
     nnflux1 = flux
;--------------------------------------------------------------------
;5) nearest_neighbors on the data itself with x, y, np
     if keyword_set(halfweights) then begin
        pixphasecorr_noisepix_halfweights, planetname, nn, apradius,chname,startaor, stopaor, /use_np 
        savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_np_hw.sav',/remove_all)
     endif else begin
        pixphasecorr_noisepix, planetname, nn, apradius,chname,startaor, stopaor, /use_np
        savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_np.sav',/remove_all)
     endelse
     restore, savename
     nnflux2 = flux
;--------------------------------------------------------------------
;6) nearest_neightbors on the data itself with x, y, xfwhm, yfwhm
     if keyword_set(halfweights) then begin
        pixphasecorr_noisepix_halfweights, planetname, nn, apradius,chname,startaor, stopaor, /use_fwhm 
        savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_fw_hw.sav',/remove_all)
     endif else begin
        pixphasecorr_noisepix, planetname, nn, apradius,chname,startaor, stopaor, /use_fwhm
        savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_fw.sav',/remove_all)
     endelse
     restore, savename
     nnflux3 = flux
     print, 'finished all nn'
  endif

;--------------------------------------------------------------------
  if keyword_set(plot_pmap) and keyword_set(pmapff) then begin
     savefile1 = strcompress(dirname + '/corrflux1_noff.sav',/remove_all)
     savefile2 = strcompress(dirname + '/corrflux2_noff.sav',/remove_all)
     savefile3 = strcompress(dirname + '/corrflux3_noff.sav',/remove_all)
     restore, savefile1 & restore, savefile2 & restore, savefile3
  endif

  if keyword_set(plot_pmap) and ~keyword_set(pmapff) then begin
     savefile1 = strcompress(dirname + '/corrflux1.sav',/remove_all)
     savefile2 = strcompress(dirname + '/corrflux2.sav',/remove_all)
     savefile3 = strcompress(dirname + '/corrflux3.sav',/remove_all)
     restore, savefile1 & restore, savefile2 & restore, savefile3
  endif


  if keyword_set(plot_nn) and ~keyword_set(halfweights) then begin
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_xy.sav',/remove_all)
     restore, savename
     nnflux1 = flux
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_np.sav',/remove_all)
     restore, savename
     nnflux2 = flux
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_fw.sav',/remove_all)
     restore, savename
     nnflux3 = flux;
  endif
  if keyword_set(plot_nn) and keyword_set(halfweights) then begin
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_xy.sav',/remove_all)
     restore, savename
     nnflux1 = flux
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_np_hw.sav',/remove_all)
     restore, savename
     nnflux2 = flux
     savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_fw_hw.sav',/remove_all)
     restore, savename
     nnflux3 = flux;
  endif

;  endelse                       ;not running the photometry
  
 
;binning
  timearr =  (planethash[aorname(startaor),'timearr'] - (planethash[aorname(startaor),'timearr'])(0))
  bmjdarr = planethash[aorname(startaor),'bmjdarr'] 
  fluxarr = planethash[aorname(startaor),'flux'] 
  numberarr = findgen(n_elements(fluxarr))
  h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  
  bin_time = dblarr(n_elements(h))
  bin_flux = bin_time
  bin_bmjd = bin_tine
  bin_corrflux1= bin_time
  bin_corrflux2= bin_time
  bin_corrflux3= bin_time
  bin_nn1 = bin_time
  bin_nn2 = bin_time
  bin_nn3 = bin_time
;  bin_corrfluxadd = bin_time
;  bin_corrfluxmult = bin_time

  c = 0
  for j = 0L, n_elements(h) - 1 do begin
     if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        bin_time[c] = median(timearr[ri[ri[j]:ri[j+1]-1]])
        bin_bmjd[c] =  median(bmjdarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, fluxarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_flux[c] = meanx     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        if keyword_set(plot_pmap) then begin
           meanclip_jk, corrflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_corrflux1[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           meanclip_jk, corrflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_corrflux2[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           meanclip_jk, corrflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_corrflux3[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        endif
        if keyword_set(plot_nn) then begin
           meanclip_jk, nnflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_nn1[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           meanclip_jk, nnflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_nn2[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           meanclip_jk, nnflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_nn3[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        endif

;        meanclip_jk, corrfluxadd[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;        bin_corrfluxadd[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;        meanclip_jk, corrfluxmult[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;        bin_corrfluxmult[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        c = c + 1
     endif
  endfor
  bin_time = bin_time[0:c-1]
  bin_bmjd = bin_bmjd[0:c-1]
  bin_flux = bin_flux[0:c-1]
  bin_corrflux1 = bin_corrflux1[0:c-1]
  bin_corrflux2 = bin_corrflux2[0:c-1]
  bin_corrflux3 = bin_corrflux3[0:c-1]
  bin_nn1 = bin_nn1[0:c-1]
  bin_nn2 = bin_nn2[0:c-1]
  bin_nn3 = bin_nn3[0:c-1]
;  bin_corrfluxadd = bin_corrfluxadd[0:c-1]
;  bin_corrfluxmult = bin_corrfluxmult[0:c-1]
  print, 'finished binning'
;--------------------------------------------------------------------
;plot the raw fluxes
;make histograms to compare widths for the different methods
;    print, 'n_elements nn', n_elements(where(finite(fluxarr))), n_elements(where(finite(nnflux1))), n_elements(where(finite(nnflux2))), n_elements(where(finite(nnflux3)))
;   print, 'n_elements corr', n_elements(where(finite(fluxarr))), n_elements(where(finite(corrflux1))), n_elements(where(finite(corrflux2))), n_elements(where(finite(corrflux3)))

;    cgPS_Open, '/Users/jkrick/irac_warm/pcrs_planets/'+planetname+'/test_pmap_hist.ps'

    cgHistoplot, fluxarr/ mean(fluxarr,/nan),/nan , thick = 2, xtitle = 'Flux', ytitle = 'Number',  $
                 datacolorname = 'black', title = planetname,/outline,xrange = [0.985, 1.015] , $
                 HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize ;, yrange = [1E-5, 1] , peak_height = 1.0
    binCenters = loc + (binsize / 2.0)
    yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
    fwhmraw = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='black', THICK=3, /OVERPLOT
    cgtext, 0.985, 3100, string(fwhmraw), color = 'black'

    if keyword_set(plot_nn) then begin

       cgHistoplot, nnflux1/ mean(nnflux1,/nan),/nan, thick = 2 , datacolorname = 'gray',/oplot,$
                    /outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
       binCenters = loc + (binsize / 2.0)
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       fwhmnn1 = 2*SQRT(2*ALOG(2))*coeff(2)
                                ;   cgPlot, binCenters, yfit, COLOR='gray', THICK=3, /OVERPLOT
       cgtext, 0.985, 2900, string(fwhmnn1), color = 'gray'
       
       cgHistoplot, nnflux2/ mean(nnflux2,/nan),/nan , thick = 2 ,datacolorname = 'red',/oplot,$
                    /outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
       binCenters = loc + (binsize / 2.0)
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       fwhmnn2 = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='red', THICK=3, /OVERPLOT
       cgtext, 0.985, 2700, string(fwhmnn2), color = 'red'
       
       cgHistoplot, nnflux3/ mean(nnflux3,/nan),/nan, thick = 2 , datacolorname = 'salmon',/oplot,$
                    /outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
       binCenters = loc + (binsize / 2.0)
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       fwhmnn3 = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='salmon', THICK=3, /OVERPLOT
       cgtext, 0.985, 2500, string(fwhmnn3), color = 'salmon'
    endif

    if keyword_set(plot_pmap) then begin
       cgHistoplot, corrflux1/ mean(corrflux1,/nan),/nan, thick = 2 , datacolorname = 'gray',/oplot,$
                    /outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
       binCenters = loc + (binsize / 2.0)
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       fwhmc1 = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='gray', THICK=3, /OVERPLOT
       cgtext, 0.985, 2300,string(fwhmc1), color = 'gray'
       
       cgHistoplot, corrflux2/ mean(corrflux2,/nan),/nan, thick = 2 ,datacolorname = 'blue',/oplot,$
                    /outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
       binCenters = loc + (binsize / 2.0)
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       fwhmc2 = 2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='blue', THICK=3, /OVERPLOT
       cgtext,0.985, 2100,string(fwhmc2), color = 'blue'
       
       cgHistoplot, corrflux3/ mean(corrflux3,/nan),/nan, thick = 2 , datacolorname = 'cyan',/oplot,$
                    /outline,HISTDATA=h, LOCATIONS=loc, BINSIZE=binsize
       binCenters = loc + (binsize / 2.0)
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       fwhmc3 =  2*SQRT(2*ALOG(2))*coeff(2)
;    cgPlot, binCenters, yfit, COLOR='cyan', THICK=3, /OVERPLOT
       cgtext, 0.985, 1900,string(fwhmc3), color = 'cyan'
    endif



;   cgPS_Close
  

    if keyword_set(plot_pmap) then begin
;pmap_correct
       p = plot(bin_time/60./60., bin_flux / mean(bin_flux), xtitle = 'Time', ytitle = 'Normalized flux', $
                '1s', title = planetname +  ' pmap correct', sym_size = 0.2,sym_filled = 1,color = 'black',$
                name = 'raw'+ string(fwhmraw))
       p1= plot(bin_time/60./60., (bin_corrflux1 / mean(bin_corrflux1,/nan)) + 0.005,  '1s', sym_size = 0.3,$
                sym_filled = 1,color = 'gray', overplot = p, name = 'X Y' + string(fwhmc1))
       p2= plot(bin_time/60./60., (bin_corrflux2 / mean(bin_corrflux2,/nan)) + 0.005,  '1s', sym_size = 0.3,$
                sym_filled = 1,color = 'blue', overplot = p, name = 'X Y NP'+ string(fwhmc2))
       p3= plot(bin_time/60./60., (bin_corrflux3 / mean(bin_corrflux3,/nan)) + 0.005,  '1s', sym_size = 0.3,$
                sym_filled = 1,color = 'cyan', overplot = p, name = 'X Y XYFWHM'+ string(fwhmc3))
       l = legend(target = [p, p1, p2, p3], position = [4, 0.999], /data, /auto_text_color)
       savefile = strcompress(dirname + 'compare_pmap_'+string(apradius)+'.png',/remove_all)
       if keyword_set(pmapff) then savefile = strcompress(dirname + 'compare_pmap_noff'+string(apradius)+'.png',/remove_all)
       p.save, savefile
    endif

;pixphasecorr_noisepix = nn
    if keyword_set(plot_nn) then begin
       p = plot(bin_time/60./60., bin_flux / mean(bin_flux), xtitle = 'Time', ytitle = 'Normalized flux', $
                title = planetname + ' nearest neighbor', '1s', sym_size = 0.2,sym_filled = 1,color = 'black', $
                name = 'raw' + string(fwhmraw))
       p4= plot(bin_time/60./60., (bin_nn1 / mean(bin_nn1,/nan)) + 0.005,  '1s', sym_size = 0.3,$
                sym_filled = 1,color = 'gray', overplot = p, name = 'X Y'+ string(fwhmnn1))
       p5= plot(bin_time/60./60., (bin_nn2 / mean(bin_nn2,/nan)) + 0.005,  '1s', sym_size = 0.3,$
                sym_filled = 1,color = 'red', overplot = p, name = 'X Y NP'+ string(fwhmnn2))
       p6= plot(bin_time/60./60., (bin_nn3 / mean(bin_nn3,/nan)) + 0.005,  '1s', sym_size = 0.3,$
                sym_filled = 1,color = 'salmon', overplot = p, name = 'X Y XYFWHM'+ string(fwhmnn3))
       l = legend(target = [p, p4, p5, p6], position = [6, 0.998], /data, /auto_text_color)
       savefile = strcompress(dirname + 'compare_nn_'+string(apradius)+'.png',/remove_all)
       if keyword_set(halfweights) then savefile = strcompress(dirname + 'compare_nn_hw'+string(apradius)+'.png',/remove_all)
       p.save, savefile
    endif





;--------------------------------------------------------------------
;make another plot with comparison between added and multiplied weights
  
;  q = plot(bin_time/60./60., bin_flux / mean(bin_flux), xtitle = 'Time', ytitle = 'Normalized flux', '1s', sym_size = 0.2,sym_filled = 1,color = 'black')
;  qa = plot(bin_time/60./60., (bin_corrfluxadd / mean(bin_corrfluxadd)) + 0.005, '1s',  sym_size = 0.2,sym_filled = 1,color = 'dark green', overplot = q)
;  qm = plot(bin_time/60./60., bin_corrfluxmult/ mean(bin_corrfluxmult) + 0.005,  '1s', sym_size = 0.2,sym_filled = 1,color = 'spring green', overplot = q)
  
  
;--------------------------------------------------------------------
;--------------------------------------------------------------------
  

    
   ;now fit these with gaussians to see what the fwhm of the gaussian model is.
    ;smaller fwhm = better
    ;or just plot normalized histograms

 end

;;look at HD7924 
;    planetname = 'HD7924b'
;    planetinfo = create_planetinfo()
;    aorname= planetinfo[planetname, 'aorname_ch2'] 
;    basedir = planetinfo[planetname, 'basedir']
;    dirname = strcompress(basedir + planetname +'/') 
;    savename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
;    restore, savename
;    nnflux1 = flux;
;    nnflux2 = flux_np
;    nnflux3 = flux_fwhm
    
;;binning
;    timearr =  (planethash[aorname(0),'timearr'] - (planethash[aorname(0),'timearr'])(0))
;    fluxarr = planethash[aorname(0),'flux'] 
;    numberarr = findgen(n_elements(nnflux1));
;    h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
    
;    bin_time = dblarr(n_elements(h))
;    bin_flux = bin_time
;    bin_corrflux1= bin_time
;    bin_corrflux2= bin_time
;    bin_corrflux3= bin_time
;    bin_nn1 = bin_time
;    bin_nn2 = bin_time
;    bin_nn3 = bin_time
;    bin_corrfluxadd = bin_time
;    bin_corrfluxmult = bin_time
    
;    c = 0
;    for j = 0L, n_elements(h) - 1 do begin
;       if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
;          bin_time[c] = median(timearr[ri[ri[j]:ri[j+1]-1]])
;          meanclip, fluxarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;          bin_flux[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;;          meanclip, corrflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;;          bin_corrflux1[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;;          meanclip, corrflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;;          bin_corrflux2[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;;          meanclip, corrflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;;          bin_corrflux3[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;          meanclip, nnflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;          bin_nn1[c] = meanx    ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;          meanclip, nnflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;          bin_nn2[c] = meanx    ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;          meanclip, nnflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;          bin_nn3[c] = meanx    ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;;          meanclip, corrfluxadd[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;;          bin_corrfluxadd[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;;          meanclip, corrfluxmult[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;;          bin_corrfluxmult[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;          c = c + 1
;       endif
;    endfor
;    bin_time = bin_time[0:c-1]
;    bin_flux = bin_flux[0:c-1]
;    bin_corrflux1 = bin_corrflux1[0:c-1]
 ;   bin_corrflux2 = bin_corrflux2[0:c-1]
 ;   bin_corrflux3 = bin_corrflux3[0:c-1]
;    bin_nn1 = bin_nn1[0:c-1]
;    bin_nn2 = bin_nn2[0:c-1]
;    bin_nn3 = bin_nn3[0:c-1]
;    bin_corrfluxadd = bin_corrfluxadd[0:c-1]
;    bin_corrfluxmult = bin_corrfluxmult[0:c-1]
;    print, 'finished binning'
;--------------------------------------------------------------------
;plot the raw fluxes
        
;pixphasecorr_noisepix = nn
;    p = plot(bin_time/60./60., bin_flux / mean(bin_flux), xtitle = 'Time', ytitle = 'Normalized flux', title = planetname + ' nearest neighbor', '1s', sym_size = 0.2,sym_filled = 1,color = 'black', name = 'raw');
;    p4= plot(bin_time/60./60., (bin_nn1 / mean(bin_nn1)) + 0.007,  '1s', sym_size = 0.3,sym_filled = 1,color = 'gray', overplot = p, name = 'X Y')
;    p5= plot(bin_time/60./60., (bin_nn2 / mean(bin_nn2)) + 0.007,  '1s', sym_size = 0.3,sym_filled = 1,color = 'red', overplot = p, name = 'X Y NP')
;    p6= plot(bin_time/60./60., (bin_nn3 / mean(bin_nn3)) + 0.007,  '1s', sym_size = 0.3,sym_filled = 1,color = 'salmon', overplot = p, name = 'X Y XYFWHM')
;    l = legend(target = [p, p4, p5, p6], position = [2.5, 0.995], /data, /auto_text_color)
    
    
    
