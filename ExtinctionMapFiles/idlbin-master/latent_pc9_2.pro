pro latent_pc9,normalize=normalize, mJy = mJy,error_bars = error_bars
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
  
  dirloc = '/Users/jkrick/irac_warm/pc9/IRAC023000/bcd/'
  aorname = ['0038163712', '0038163968','0038163200', '0038163456', '0038162688', '0038162944', '0038164224', '0038164480','0038166272','0038166528','0038165248','0038165504','0038164736','0038164992','0038165760','0038166016']
;  xcen =[88.04,162.0,42.,201.0,120.0,143.0,59.,218.]
;  ycen = [133.69,115.48,116.,143.0,202., 31.0,160.,188.]
  xcen =[88.04,162.0,42.,201.0,120.0,59.,218.]
  ycen = [133.69,115.48,116.,143.0,202.,160.,188.]
;  roughmag = [3, 4, 1,2,0,4,1,2]
;  colorarr  = [bluecolor, purplecolor, yellowcolor, greencolor, redcolor, purplecolor, yellowcolor, greencolor]
  roughmag = [3, 4, 1,2,0,1,2]
  colorarr  = [bluecolor, purplecolor, yellowcolor, greencolor, redcolor, yellowcolor, greencolor]
  

;for ch2
  
;  aorname = ['0038166272','0038166528','0038165248','0038165504','0038164736','0038164992'];,'0038165760','0038166016']  


;because idl defines from 0,0 not 1,1
  xcen = xcen - 1
  ycen = ycen - 1
  
;set up for plotting everything together
;the standard is to leave units in electrons
  yr = [1E2,1E4]
  yt = 'Electrons'
  
  if keyword_set(NORMALIZE) eq 1 then begin
     yr = [1E-5, 1]
     yt = 'Normalized Flux'
  endif 
  
  if keyword_set(mJy) eq 1 then begin
     yr = [1E0,1E3]       
     yt = 'mJy'
  endif
  
  plot, findgen(10), findgen(10), yrange = yr, xrange = [0,1], /nodata, title = 'Ch 1', /ylog,  ytitle = yt, xtitle = 'Time in hours', ystyle = 1, xstyle = 1
  
;now start amassing photometry
  
  for n = 0, 1 do begin; n_elements(xcen) - 1 do begin ; for each star location
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
        print, 'working on air', aorname(a)

        CD,dirloc+string(aorname(a))
        
        command1 =  ' ls IRAC.1*bcd_fp.fits >  bcd.txt'
        spawn, command1
        command1 =  ' ls IRAC.1*bcd_fp_SatFixed.fits >>  bcd.txt'
        spawn, command1
        
        readcol,strcompress(dirloc+string(aorname(a))+ '/bcd.txt', /remove_all), bcdlist, format='A', /silent 
        
        for i = 1, n_elements(bcdlist) - 1  do begin ;ignore the first frame

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
           xbkg_radius = 40 < (254 - xcen(n))
           ybkg_radius = 40 < (254 - ycen(n))

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
           
           
;fit a gaussian to this distribution
           binsize = 5.
           plothist, bkgd, xhist, yhist, bin = binsize,/noplot
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

           if keyword_set(mJy) eq 1 then begin
                                ;convert back to Mjy/sr
              flux = flux /sbtoe
              fluxerr = fluxerr /sbtoe
              sky = sky / sbtoe
              sigma = sigma/sbtoe
                                ;convert into microjy
              flux = flux * 34.9837
              fluxerr = fluxerr * 34.9837 
              sky = sky * 34.9837 
              sigma  = sigma * 34.9837
           endif
           
           flux1arr(count) = flux - result(0)
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
     print, 'start times', timearr(0)
     print, 'flux1arr', flux1arr

     if keyword_set(NORMALIZE) eq 1 then begin
        flux1err_percent = flux1err / flux1arr
        sky_percent = skyarr / flux1arr
        
        flux1arr = flux1arr / max(flux1arr)
        flux1err = flux1err_percent*flux1arr
        
        skyarr = sky_percent*flux1arr
     endif 
     
  
     if keyword_set(error_bars) eq 1 then begin
        oploterror, timearr, flux1arr,flux1err, color = colorarr(n), errcolor = colorarr(n)
        oplot, timearr, skyarr, color = colorarr(n), linestyle = 2
     endif else begin
        oplot, timearr, flux1arr, color = colorarr(n)
        oplot, timearr, skyarr , color = colorarr(n),linestyle = 2
;        oplot, timearr, flux1arr,  color = colorarr(n), thick = 3
     endelse
     
     

     
  endfor                        ;for each star
  
;ps_close, /noprint,/noid
  
end
