pro plot_noisepix, planetname, bin_level

;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  chname = planetinfo[planetname, 'chname']
  ra_ref = planetinfo[planetname, 'ra']
  dec_ref = planetinfo[planetname, 'dec']
  aorname = planetinfo[planetname, 'aorname']
  basedir = planetinfo[planetname, 'basedir']
  utmjd_center =  planetinfo[planetname, 'utmjd_center']
  transit_duration =  planetinfo[planetname, 'transit_duration']
  period =  planetinfo[planetname, 'period']
  intended_phase = planetinfo[planetname, 'intended_phase']
  stareaor = planetinfo[planetname, 'stareaor']
  plot_norm= planetinfo[planetname, 'plot_norm']
  plot_corrnorm = planetinfo[planetname, 'plot_corrnorm']
  
  dirname = strcompress(basedir + planetname +'/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'.sav')
  restore, savefilename
  
;-------------------------------------------------------
;make a bunch of plots before binning
;-------------------------------------------------------


;-----
  for a = 0 , n_elements(aorname) -1 , 2 do begin                        
     print, 'working on aorname ', aorname(a)
  
;  am = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'xcen'],'1-',  ytitle = 'X pix', title = planetname, xtitle = 'Time(hrs)', xrange = [2.5,3.5], yrange = [15.0, 15.2]) ;
;  am.save, dirname +'x_time_ch'+chname+'_part.png'
  
;------
  
;  ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'ycen'],'1-', ytitle = 'Y pix', xtitle = 'Time(hrs)', xrange = [2.5,3.5], title = planetname, yrange = [14.8, 15.0])
;  ay.save, dirname +'y_time_ch'+chname+'_part.png'
  
;------
; ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'flux'],'1-',ytitle = 'Flux', xtitle = 'Time(hrs)', xrange = [2.5,3.5], title = planetname, yrange = [1.42, 1.5])
;  ay.save, dirname +'raw_flux_time_ch'+chname+'_part.png'
  
;------
;  ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'bkgd'],'1-',  ytitle = 'Bkgd', xtitle = 'Time(hrs)', xrange = [2.5,3.5], title = planetname)
;  ay.save, dirname +'bkgd_time_ch'+chname+'_part.png'
  
;------
;  ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'np'],'1-', ytitle = 'NP', xtitle = 'Time(hrs)', xrange = [2.5,3.5], title = planetname, yrange = [4.0, 5.0])
;  ay.save, dirname +'np_time_ch'+chname+'_part.png'
  
;------
  
;can I see sigma position changing with time on timescales of 0.05hours = 3 min? 
;at 8.5s from one frame start to the next, that is 1330 individual frames (21 cubes) 
;specifically at the time of the spike in noise pix?

;  bin_level = 1330L
;

; binning
     numberarr = findgen(n_elements((planethash[aorname(a),'timearr'])))
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     print, 'omin', om, 'nh', n_elements(h)
  
  
;mean together the flux values in each phase bin     
     bin_xcen = dblarr(n_elements(h))
     bin_ycen = dblarr(n_elements(h))
     bin_xsigma = dblarr(n_elements(h))
     bin_ysigma = dblarr(n_elements(h))
     bin_bkgd = dblarr(n_elements(h))
     bin_flux = dblarr(n_elements(h))
     bin_timearr = dblarr(n_elements(h))
     bin_np = dblarr(n_elements(h))

     timearr = (planethash[aorname(a),'timearr'] - (planethash[aorname(a),'timearr'])(0))/60./60.
     
     
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
        
;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
           
           meanclip, (planethash[aorname(a),'xcen'])[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_xcen[c] = meanx  ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           bin_xsigma[c]=sigmax
           
           meanclip,  (planethash[aorname(a),'ycen'])[ri[ri[j]:ri[j+1]-1]], meany, sigmay
           bin_ycen[c] = meany  ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           bin_ysigma[c] = sigmay
           
           meanclip,  (planethash[aorname(a),'bkgd'])[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
           bin_bkgd[c] = meansky ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip,  (planethash[aorname(a),'flux'])[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip,  (planethash[aorname(a),'np'])[ri[ri[j]:ri[j+1]-1]], meannp, sigmanp
           bin_np[c] = meannp   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip,  timearr[ri[ri[j]:ri[j+1]-1]], meantimearr, sigmatimearr
           bin_timearr[c]=meantimearr
           
           c = c + 1
        endif
     endfor
     
     bin_xcen = bin_xcen[0:c-1]
     bin_ycen = bin_ycen[0:c-1]
     bin_xsigma = bin_xsigma[0:c-1]
     bin_ysigma = bin_ysigma[0:c-1]
     bin_bkgd = bin_bkgd[0:c-1]
     bin_flux = bin_flux[0:c-1]
     bin_timearr = bin_timearr[0:c-1]
     bin_np = bin_np[0:c-1] 
     
;print, 'bin_xcen', bin_xcen;
     print, 'bin)timearr', min(bin_timearr), max(bin_timearr)

;------
;plot the binned stuff
     setxrange = [0.,8.]
     
     pp = plot(bin_timearr, bin_xsigma, '1-', title = planetname, $
               xtitle = 'Time(hrs)', ytitle = 'X sigma', $;sym_size = 0.2, sym_filled =1,$
               xrange = setxrange, yrange = [0,0.03])
     pp.save, dirname +'xsigma_time_'+aorname(a)+'.png'

     pq = plot(bin_timearr, bin_ysigma, '1-', title = planetname,  $
               xtitle = 'Time(hrs)', ytitle = 'Y sigma', $;sym_size = 0.2, sym_filled =1,$
               xrange = setxrange)
     pq.save, dirname +'ysigma_time_'+aorname(a)+'.png'

;     pr = plot(bin_timearr, bin_flux/plot_norm, '1-',  $
;               xtitle = 'Time(hrs)',$;sym_size = 0.2, sym_filled =1, $
;               ytitle = 'Normalized Flux', title = planetname, yrange = [0.995,1.005], $
;               xrange = setxrange) 
;     pr.save, dirname +'flux_time_'+aorname(a)+'.png'

     ps = plot(bin_timearr, bin_np, '1-', title = planetname,  $
               xtitle = 'Time(hrs)', ytitle = 'NP', $;sym_size = 0.2, sym_filled =1,$
               xrange = setxrange)
     ps.save, dirname +'np_time_'+aorname(a)+'.png'

    ps = plot(bin_timearr, bin_xcen, '1-', title = planetname,  $
               xtitle = 'Time(hrs)', ytitle = 'xcen', $;sym_size = 0.2, sym_filled =1,$
               xrange = setxrange)
     ps.save, dirname +'xcen_time_'+aorname(a)+'.png'

    ps = plot(bin_timearr, bin_ycen, '1-', title = planetname,  $
               xtitle = 'Time(hrs)', ytitle = 'ycen', $;sym_size = 0.2, sym_filled =1,$
               xrange = setxrange)
     ps.save, dirname +'ycen_time_'+aorname(a)+'.png'


;-------------
;second order
;sigma sigma

;re- binning
     numberarr = findgen(n_elements(bin_xcen))
     h = histogram(numberarr, OMIN=om, binsize = 10L, reverse_indices = ri)
     print, 'omin', om, 'nh', n_elements(h)
  
  
;mean together the flux values in each phase bin     
     bin_xsigmasigma = dblarr(n_elements(h))
     bin_ysigmasigma = dblarr(n_elements(h))
     bin_sigmatimearr = dblarr(n_elements(h))

     c = 0
     for j = 0L, n_elements(h) - 1 do begin
        
;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
           
           meanclip, (bin_xsigma)[ri[ri[j]:ri[j+1]-1]], meanxsigma, sigmasigmax
           bin_xsigmasigma[c]=sigmasigmax
           
           meanclip, (bin_ysigma)[ri[ri[j]:ri[j+1]-1]], meanysigma, sigmasigmay
           bin_ysigmasigma[c]=sigmasigmay

           meanclip,  bin_timearr[ri[ri[j]:ri[j+1]-1]], meantimearr, sigmatimearr
           bin_sigmatimearr[c]=meantimearr
           
           c = c + 1
        endif
     endfor
     
     bin_xsigmasigma = bin_xsigmasigma[0:c-1]
     bin_ysigmasigma = bin_ysigmasigma[0:c-1]
     bin_sigmatimearr = bin_sigmatimearr[0:c-1]
;---------
;plots

;     pp = plot(bin_sigmatimearr, bin_xsigmasigma, '1-', title = planetname, $
;               xtitle = 'Time(hrs)', ytitle = 'X sigma sigma', $;sym_size = 0.2, sym_filled =1,$
;               xrange = setxrange)
;     pp.save, dirname +'xsigmasigma_time_'+aorname(a)+'.png'

;     pq = plot(bin_sigmatimearr, bin_ysigmasigma, '1-', title = planetname,  $
;               xtitle = 'Time(hrs)', ytitle = 'Y sigma sigma', $;sym_size = 0.2, sym_filled =1,$
;               xrange = setxrange)
;     pp.save, dirname +'ysigmasigma_time_'+aorname(a)+'.png'


  endfor  ; for each AOR

  
END
  
  
