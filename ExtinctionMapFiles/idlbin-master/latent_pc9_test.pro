pro latent_pc9_test
  !P.multi =[0,1,1]
;  ps_open, filename='/Users/jkrick/iwic/latent/above_back.ps',/portrait,/square,/color


  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
  whitecolor = FSC_COLOR("white", !D.Table_Size-9)
  colorarr = [ redcolor ,bluecolor,greencolor,yellowcolor,cyancolor,orangecolor,purplecolor,whitecolor]
  
 
;for ch1
  
  dirloc = '/Users/jkrick/irac_warm/latent/'
  aorname ='r32898816';'r38249984'; 'r38348032'
  xcen = 191.6;48.8;69.3
  ycen = 186.0; 6.8; 80.2

;because idl defines from 0,0 not 1,1
  xcen = xcen - 1
  ycen = ycen - 1
  
;set up for plotting everything together
;the standard is to leave units in electrons
  yr = [1E2,1E4]
  yt = 'Electrons'
  
  plot, findgen(10), findgen(10), yrange = yr, xrange = [0,3], /nodata, title = 'Ch 1', /ylog,  ytitle = yt, xtitle = 'Time in hours', ystyle = 1, xstyle = 1
  
;now start amassing photometry
  
  flux1arr =fltarr(165*10)
  flux1err =fltarr(165*10)
  skyarr = fltarr(165*10)
  timearr = dblarr(165*10)
  time = double(0)
  
  count = 0
  
  CD,dirloc+aorname + '/ch1/bcd/'
  command1 =  ' ls SPITZER_I1*_bcd.fits >  bcd.txt'
  spawn, command1
  command1 =  ' ls SPITZER_I1*_bcd_SatFixed.fits >>  bcd.txt'
  spawn, command1
        
  readcol, dirloc+aorname+ '/ch1/bcd/bcd.txt', bcdlist, format='A', /silent 
        
  for i = 91, n_elements(bcdlist) - 1  do begin ;ignore the first frame   first 2 frames if using the saturation corrected value.
     print, 'working on bcdlist', bcdlist(i), n_elements(bcdlist)
     fits_read, bcdlist(i), data, imheader
                                ;get rid of NaN's
     data[where(finite(data) lt 1)] = 0
           
                                ;convert the data from MJy/sr to electrons
     gain = 3.7
     exptime =  sxpar(imheader, 'EXPTIME') 
     fluxconv =  sxpar(imheader, 'FLUXCONV') 
     sbtoe = gain*exptime/fluxconv
     data = data*sbtoe
           
                                ;instead try doing this manually
     sq2 =2    
     flux = total(data[xcen - sq2:xcen+sq2, ycen-sq2:ycen + sq2])
     fluxerr = stddev(data[xcen - sq2:xcen+sq2, ycen-sq2:ycen + sq2])

;randomly choose centers to look at background regions
     nrand = 100
     xbkg_radius = 36 < (xcen -1);(254 - xcen)
     ybkg_radius = 36 < (ycen - 1);(254 - ycen)

     x = randomu(S,nrand) *(xbkg_radius * 2) + xcen - xbkg_radius
     y = randomu(S,nrand) *(ybkg_radius * 2) + ycen - ybkg_radius
     bkgd = fltarr(nrand)
     for ti = 0, n_elements(x) - 1 do begin
        pixval = data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2]
        ga = where(finite(pixval) gt 0)
        if n_elements(ga) gt .85*n_elements(pixval) then  bkgd(ti) =  total(data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2],/nan) ; / n_elements(ga)
     endfor
     ta = where(bkgd ne 0)
     bkgd = bkgd(ta)
           
;fit a gaussian to this distribution
     binsize = 200              ; 5.
     plothist, bkgd, xhist, yhist, bin = binsize,/noplot
     start = [1000.,200., 2.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                                       ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet)       ;./quiet    

;now plot the histogram of stddevs of all unmasked pixels
;           plothist, bkgd, xjunk, yjunk, bin = binsize, xtitle = 'Background (electrons)', ytitle='Number', $
;                     xrange=[-8*result(1) + result(0),result(0) + 8*result(1)] 
;plot the fitted gaussian and print results to plot
;           xarr = findgen(1000)/ 10.          
;           oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
;                  thick=5, linestyle =  2 ;color = greencolor,
                     
                                ;keep track of time of observation
     time= fxpar(imheader, 'SCLK_OBS')
;           print, 'flux in electrons', flux, sky, flux-sky

           
     flux1arr(count) = flux - result(0)
     print, 'flux', flux - result(0)
     flux1err(count) = fluxerr
     skyarr(count) = 2.*result(1)
     timearr(count) = time
     count = count + 1
           
  endfor                        ;For Each image in the bcdlist
        
     
     flux1arr = flux1arr[0:count - 1]
     flux1err = flux1err[0:count - 1]
     skyarr = skyarr[0:count-1]
     timearr= timearr[0:count - 1]
     
     timearr = (timearr - timearr(0))/3600.

     print, 'pct down for star ',  flux1arr(0) / flux1arr(1), flux1arr(0), flux1arr(1)
  
                                ;binning in time
     h = histogram(timearr, reverse_indices=ri, omin=om,binsize = 0.1) ;0.08
     flux1bin = fltarr(n_elements(h))
     timebin = findgen(n_elements(h))*0.1
     err1bin = fltarr(n_elements(h))
     skybin = fltarr(n_elements(h))
     
     for j = 0, n_elements(h) - 1 do begin
        meanclip,flux1arr[ri[ri[j]:ri[j+1]-1]], meanbin, sigmabin, clipsig = 2
        flux1bin[j] =  meanbin
        err1bin = sigmabin
        
        meanclip,skyarr[ri[ri[j]:ri[j+1]-1]], skymean, skysigma, clipsig = 2
        skybin[j] =  skymean
        
     endfor                     ;end binning
     
     if keyword_set(error_bars) eq 1 then begin
        oploterror, timebin, flux1bin,flux1err, color = colorarr(n), errcolor = colorarr(n)
        oplot, timebin, skybin, color = colorarr(n), linestyle = 2
     endif else begin
        oplot, timebin, flux1bin, color = colorarr(0)
        oplot, timebin, skybin , color = colorarr(0),linestyle = 2
;        oplot, timearr, flux1arr,  color = colorarr(n), thick = 3
     endelse
     
     
;try fitting
;try just after 20 min
     a = where(timebin gt 0.15)
     timebin = timebin(a)
     flux1bin = flux1bin(a)   
     
     start = [1000.,1.]
     noise = fltarr(n_elements(flux1bin))
     noise[*] = 1                                                   ;equally weight the values
     result= MPFITFUN('exponential',timebin,flux1bin, noise, start,/quiet) ;./quiet    
     print, 'results', result(0), result(1)
;     oplot, findgen(1000)/10, result(0)*exp(-(findgen(1000)/10)/(result(1))), linestyle = 3, color = colorarr(n)
     oplot, findgen(1000)/10, result(0)*exp(-(findgen(1000)/10)/(4.4)), linestyle = 3
     
;ok, now when will this reach 2xnoise in the background?
     
                                ;assume we pull a fits file at time = 4hrs
     latent_flux = result(0)*exp(-4/result(1))
;     print, 'latent flux', latent_flux
     
;now ask where does the latent flux eq the sky noise background.
     edge = where(flux1bin le skybin, edgecount)
     if edgecount gt 1 then print, 'times', timebin(edge)
     
;and ask where the fit crosses the noise
     fit_edge = where( result(0)*exp(-(timebin)/(4.4)) le skybin,fitcount)
     if fitcount gt 1 then print, 'times for fit', timebin(fit_edge)

 
end
