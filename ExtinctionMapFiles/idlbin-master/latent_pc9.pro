pro latent_pc9,normalize=normalize, mJy = mJy,error_bars = error_bars
  !P.multi =[0,1,1]
  ;ps_open, filename='/Users/jkrick/iwic/latent/above_back.ps',/portrait,/square,/color
ps_start, filename= '/Users/jkrick/iwic/latent/above_back_ch1.ps'




  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
  whitecolor = FSC_COLOR("white", !D.Table_Size-9)
  colorarr = [ redcolor ,bluecolor,greencolor,orangecolor,cyancolor,orangecolor,purplecolor,whitecolor]
  
  ra = [86.983887,86.990761, 82.220307, 82.215477]
  dec = [24.686835, 24.661534, 32.477487, 32.422855]

;for ch1
  
  dirloc = '/Users/jkrick/irac_warm/pc9/IRAC023000/bcd/'
  aorname = ['0038163712', '0038163968','0038163200', '0038163456', '0038162688', '0038162944', '0038164224', '0038164480','0038166272','0038166528','0038165248','0038165504','0038164736','0038164992','0038165760','0038166016']
;  xcen =[88.04,162.0,42.,201.0,120.0,143.0,59.,218.]
;  ycen = [133.69,115.48,116.,143.0,202., 31.0,160.,188.]
  xcen =[88.04,162.0,42.,201.0,120.0,59.,218.]
  ycen = [133.69,115.48,116.,143.0,202.,160.,188.]
;  roughmag = [3, 4, 1,2,0,4,1,2]
;  colorarr  = [bluecolor, purplecolor, yellowcolor, greencolor, redcolor, purplecolor, yellowcolor, greencolor]
  roughmag = [3, 4, 1,2,0,1,2]
  colorarr  = [bluecolor, purplecolor, orangecolor, greencolor, redcolor, orangecolor, greencolor]
  

;for ch2
  
;  aorname = ['0038166272','0038166528','0038165248','0038165504','0038164736','0038164992'];,'0038165760','0038166016']
;  xcen =[83.76,158.7,38.0,198.1,120.0,59.,218.]
;  ycen = [133.69,115.2,115.6,142.7,202.,160.,188.]

;because idl defines from 0,0 not 1,1
  xcen = xcen - 1
  ycen = ycen - 1
  
;set up for plotting everything together
;the standard is to leave units in electrons
  yr = [1E2,1E4]
;  yr = [-400, 400]
  yt = 'Electrons'
  
  if keyword_set(NORMALIZE) eq 1 then begin
     yr = [1E-5, 1]
     yt = 'Normalized Flux'
  endif 
  
  if keyword_set(mJy) eq 1 then begin
     yr = [1E0,1E3]       
     yt = 'mJy'
  endif
  
  ;plot, findgen(10), findgen(10), yrange = y, xrange = [0,1], /nodata, title = 'Ch 2',  ytitle = yt, xtitle = 'Time in hours', ystyle = 1, xstyle = 1
  plot, findgen(10), findgen(10), yrange = yr, xrange = [0,8], /nodata, title = 'Ch 1 Latent Profile',  ytitle = yt, xtitle = 'Time in hours', ystyle = 1, xstyle = 1,/ylog

  legend, ['1st mag. star', 'Background', 'Model'], linestyle = [0,2,3], position = [5,8000]
;now start amassing photometry
  
  for n = 2, 2 do begin; n_elements(xcen) - 1 do begin ; for each star location
     print, 'working on star', n
     flux1arr =fltarr(165*10)
     flux1err =fltarr(165*10)
     skyarr = fltarr(165*10)
     timearr = dblarr(165*10)
     time = double(0)
     
     count = 0
     
     case n of 
        0: s=0
        1: s=0
        2: s=2
        3: s=2
        4: s=4
        5: s=4
        6: s=6
        7: s=6
     endcase
     
     for a = s, n_elements(aorname) - 1 do begin
       print, 'working on', aorname(a), n_elements(aorname)
        CD,dirloc+string(aorname(a))
        
        command1 =  ' ls IRAC.1*bcd_fp.fits >  bcd.txt'
        spawn, command1
        command1 =  ' ls IRAC.1*bcd_fp_SatFixed.fits >>  bcd.txt'
        spawn, command1
        
        readcol,strcompress(dirloc+string(aorname(a))+ '/bcd.txt', /remove_all), bcdlist, format='A', /silent 
        
        for i = 2, n_elements(bcdlist) - 1  do begin ;ignore the first frame   first 2 frames if using the saturation corrected value.
           print, 'working on bcdlist', bcdlist(i)
           fits_read, bcdlist(i), data, imheader
                                ;get rid of NaN's
           data[where(finite(data) lt 1)] = 0
           
                                ;convert the data from MJy/sr to electrons
           gain = 3.7
           exptime =  sxpar(imheader, 'EXPTIME') 
           fluxconv =  sxpar(imheader, 'FLUXCONV') 
           sbtoe = gain*exptime/fluxconv
           data = data*sbtoe
           
                                ;do photometry
;           aper, data, xcen(n), ycen(n), flux, fluxerr, sky, skyerr, 1, [3], [7,10],/NAN,/flux,/silent,/exact
;           print, 'aper flux', flux, sky*!Pi*2^2, format = '(A, F10.2,F10.2)'
           
                                ;instead try doing this manually
           sq2 =2    
           flux = total(data[xcen(n) - sq2:xcen(n)+sq2, ycen(n)-sq2:ycen(n) + sq2])
           fluxerr = stddev(data[xcen(n) - sq2:xcen(n)+sq2, ycen(n)-sq2:ycen(n) + sq2])
;           sky1 = total(data[xcen(n) +12:xcen(n) + 16, ycen(n) +8:ycen(n) + 1sq2])
           
;want to know what the noise in the local background is.
;take 1
;           skydata = [data[xcen(n) +12:xcen(n) + 20, ycen(n) -20:ycen(n) + 20] , data[xcen(n) -20:xcen(n) - 12, ycen(n) -20:ycen(n) + 20]];
;,  data[xcen(n) -12:xcen(n) + 12, ycen(n) +12:ycen(n) +20] ];,  data[xcen(n) -12:xcen(n) + 12, ycen(n)-12:ycen(n) -20]  ]
;           mmm, skydata, sky, sigma, skew
;           sky = sky *25.       ; number of pixels in aperture above
           
;take 2
;randomly choose centers to look at background regions
           nrand = 100
           xbkg_radius = 36 < (254 - xcen(n))
           ybkg_radius = 36 < (254 - ycen(n))

           x = randomu(S,nrand) *(xbkg_radius * 2) + xcen(n) - xbkg_radius
           y = randomu(S,nrand) *(ybkg_radius * 2) + ycen(n) - ybkg_radius
           bkgd = fltarr(nrand)
           for ti = 0, n_elements(x) - 1 do begin
              pixval = data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2]
              ga = where(finite(pixval) gt 0)
              if n_elements(ga) gt .85*n_elements(pixval) then  bkgd(ti) =  total(data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2],/nan); / n_elements(ga)
           endfor
           ta = where(bkgd ne 0)
           bkgd = bkgd(ta)
           
           if keyword_set(mJy) eq 1 then begin
                                ;convert back to Mjy/sr
              flux = flux /sbtoe
              fluxerr = fluxerr /sbtoe
;              sky = sky / sbtoe
;              sigma = sigma/sbtoe
              bkgd = bkgd / sbtoe
                                ;convert into microjy
              flux = flux * 34.9837
              fluxerr = fluxerr * 34.9837 
;              sky = sky * 34.9837 
;              sigma  = sigma * 34.9837
              bkgd = bkgd * 34.9837
           endif
;;           print, 'bkgd', bkgd
;fit a gaussian to this distribution
           binsize = 100; 5.
;           print, 'starting plothist'
           plothist, bkgd, xhist, yhist, bin = binsize,/noplot
;           print, 'finished plothist'
;           print,'flux, mean, stddev , simple sky, sigma', flux, mean(bkgd), stddev(bkgd), sky, sigma
           start = [1000.,200., 2.]
           noise = fltarr(n_elements(yhist))
           noise[*] = 1                                                 ;equally weight the values
           result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ;./quiet    

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

           print, 'flux', flux, result(0)
           flux1arr(count) = flux - result(0)
;           print, 'flux', flux - result(0)
           flux1err(count) = fluxerr
           skyarr(count) = 2.*result(1)
           timearr(count) = time
           count = count + 1
           
        endfor                  ;For Each image in the bcdlist
        
     endfor                     ;for each aorname
     
     flux1arr = flux1arr[0:count - 1]
     flux1err = flux1err[0:count - 1]
     skyarr = skyarr[0:count-1]
     timearr= timearr[0:count - 1]
     
     timearr = (timearr - timearr(0))/3600.
;     print, 'start times', timearr(0)
     print, 'flux1arr', flux1arr


     print, 'pct down for star ', roughmag(n) , flux1arr(0) / flux1arr(1), flux1arr(0), flux1arr(1)

     if keyword_set(NORMALIZE) eq 1 then begin
        flux1err_percent = flux1err / flux1arr
        sky_percent = skyarr / flux1arr
        
        flux1arr = flux1arr / max(flux1arr)
        flux1err = flux1err_percent*flux1arr
        
        skyarr = sky_percent*flux1arr
     endif 
     
                                ;binning in time
     h = histogram(timearr, reverse_indices=ri, omin=om,binsize = 0.08) ;0.05
     flux1bin = fltarr(n_elements(h))
     timebin = findgen(n_elements(h))*0.08
     err1bin = fltarr(n_elements(h))
     skybin = fltarr(n_elements(h))
     
     for j = 0, n_elements(h) - 1 do begin
        meanclip,flux1arr[ri[ri[j]:ri[j+1]-1]], meanbin, sigmabin, clipsig = 2
        flux1bin[j] =  meanbin
        err1bin = sigmabin
        
        meanclip,skyarr[ri[ri[j]:ri[j+1]-1]], skymean, skysigma, clipsig = 2
        skybin[j] =  skymean
        
     endfor                     ;end binning
     
;     print, 'flux1bin', flux1bin
;     print, 'skybin', skybin
     if keyword_set(error_bars) eq 1 then begin
        oploterror, timebin, flux1bin,flux1err;, color = colorarr(n), errcolor = colorarr(n)
        oplot, timebin, skybin, linestyle = 2;, color = colorarr(n)
     endif else begin
        oplot, timebin, flux1bin, thick = 3;, color = colorarr(n)
        oplot, timebin, skybin ,linestyle = 2;, color = colorarr(n)
        oplot, timebin, -skybin ,linestyle = 2;, color = colorarr(n)
;        oplot, timearr, flux1arr, color = colorarr(n)
;        oplot, timearr, skyarr , color = colorarr(n),linestyle = 2
;        oplot, timearr, -skyarr , color = colorarr(n),linestyle = 2
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
     ;oplot, findgen(1000)/10, result(0)*exp(-(findgen(1000)/10)/(result(1))), linestyle = 3, color = colorarr(n)
 
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

;     save, filename = '/Users/jkrick/irac_warm/mag3_fluxbin.sav', flux1bin
;     save, filename = '/Users/jkrick/irac_warm/mag3_skybin.sav', skybin
;     save, filename = '/Users/jkrick/irac_warm/mag3_timebin.sav', timebin


     
  endfor                        ;for each star
  
;ps_close, /noprint,/noid
  ps_end, /png

end
