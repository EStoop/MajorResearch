pro size_sigma
nan = alog10(-1)
ch = 0

ps_start, filename= '/Users/jkrick/Virgo/IRAC/test_size_sigma_ch2.ps'
!P.multi = [0,1,1]

if ch eq 0 then begin
;   fits_read, '/Users/jkrick/Virgo/IRAC/ch1_drizzle/ch1_Combine-mosaic/mosaic_bkgd.fits', data, header
  fits_read, '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/mosaic_bkgd.fits', data, header
   ;fits_read, '/Users/jkrick/Virgo/IRAC/s18p14/ch1/ch1_Combine-mosaic/mosaic_cov.fits', covdata, covheader
   cd, '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic'
   command = 'sex mosaic_bkgd.fits -c ../../default.sex'
  ; spawn, command
   fits_read, '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/segmentation.fits', segdata, segheader
   print, '****ch1*****'

endif

if ch eq 1 then begin
   fits_read, '/Users/jkrick/Virgo/IRAC/ch2/ch2_Combine-mosaic/mosaic_bkgd.fits', data, header
   ;fits_read, '/Users/jkrick/Virgo/IRAC/s18p14/ch2/ch2_Combine-mosaic/mosaic_cov.fits', covdata, covheader
   cd, '/Users/jkrick/Virgo/IRAC/ch2/ch2_Combine-mosaic'
;   command = 'sex mosaic_bkgd.fits -c ../../default.sex'
;   spawn, command
   fits_read, '/Users/jkrick/Virgo//IRAC/ch2/ch2_Combine-mosaic/segmentation.fits', segdata, segheader
   print, '****ch2*****'
endif

;try converting the data image into electrons.
;sbtoe = gain*exptime/flux_conv
;Mjypersr * sbtoe = electrons 

;sbtoe = (covdata*100) * 3.7 / .1198
;data = data * sbtoe

;first use the segmentation image as a mask
   segmask = where(segdata gt 0)
   data(segmask) = nan
;;   fits_write, strcompress('/Users/jkrick/Virgo/IRAC/s18p14/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_seg.fits',/remove_all), data, header
print, 'how many masked pixels are there', n_elements(segmask), n_elements(data)

if ch eq 0 then begin
xcen = [996,668,1427,1338,591,1391,2369,3363,1149]
ycen = [5399,4427,4323,3560,3034,2938,1852,1836,949]
endif
if ch eq 1 then begin
xcen = [716,664,755,755,1687,3128];4805
ycen = [2952,3574,4337,937,1894,2423];3711
endif 
e = 1000
sigmaarr = fltarr(e)
sizearr = fltarr(e)

for n = 1, e do begin 
   size = 1*n
   sumarr = fltarr(n_elements(xcen))
   for i = 0, n_elements(xcen) - 1 do begin 
      goodarea = where(finite(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.]) gt 0, goodcount)
      sumarr(i) = total(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.],/nan) / goodcount
   endfor
   
   sigmaarr(n-1) = stddev(sumarr)
   sizearr(n-1) = (sqrt(goodcount))*0.6
endfor

plot, sizearr, sigmaarr, thick = 3, charthick = 3,/ylog, ytitle = 'Noise in Mjy/sr/pixel', xtitle = 'Binning Length (arcsec on bottom, native pixels on top)',  psym = 2, /xlog, xstyle = 9, xrange = [1,1000], xthick = 3, ythick =3, ystyle = 9, yrange = [.0001, .01]

axis, xaxis = 1, xrange = [1/1.2,1000/1.2], xstyle = 1,/xlog, xthick = 3
if ch eq 0 then axis, yaxis =1, ystyle = 1, /ylog, ythick = 3, yrange = [.0001*1E6*1E-26*8.33E13*1E9, .01*1E6*1E-26*8.33E13*1E9]
if ch eq 1 then axis, yaxis =1, ystyle = 1, /ylog, ythick = 3, yrange = [.0001*1E6*1E-26*6.67E13*1E9, .01*1E6*1E-26*6.67E13*1E9]

;try fitting a function
;start = [0.009]
;noise = fltarr(n_elements(sigmaarr))
;noise[*] = 1.
;result= MPFITFUN('power',sizearr,sigmaarr, noise, start)    ;ICL
;oplot, sizearr,  result(0)*(sizearr)^(-0.5) 

oplot, sizearr, 0.0035/ sizearr, thick = 3


;----------------------------------------------------
;repeat with a larger masked version of the image

;first make the larger masks
;try increasing the mask size of all objects uniformly
   readcol, 'test.cat', NUMBER                 , X_WORLD                , Y_WORLD                , X                , Y                , MAG_AUTO               , FLUX_BEST              , FLUXERR_BEST           , MAG_BEST               , MAGERR_BEST            , BACKGROUND             , FLUX_MAX               , ISOAREA_IMAGE          , ALPHA_J2000            , DELTA_J2000            , A_IMAGE                , B_IMAGE                , THETA_IMAGE            , MU_THRESHOLD           , MU_MAX                 , FLAGS                  , FWHM_IMAGE             , CLASS_STAR             , ELLIPTICITY            

   ;print, ellipticity(24)

   sma = sqrt(isoarea_image/(!PI*(1-ellipticity)))
   smb = isoarea_image/(!PI*sma)
   sma = 1.5*sma                ; 3.0,2.3 ,1.6                     3.2;3.0*sma;2.6*sma;2*sma
   smb = 1.5*smb                ; 3.0,2.3,1.6                     3.2;3.0*smb;2.6*smb;2*smb
   
;convert theta from degres into radians
   theta = theta_image*(!PI/180)
   increment =  (2*!PI)/500
  
   xmax = 6261                  ;1601;2049;1900;2866;1900;2866;2048;2000;2800;1000
   ymax = 6155                  ;1801;3148;2100;3486;2100;3486;3148;2000;2700;1000
   
   ;for each detection
   for d = long(0), n_elements(number) - 1 do begin

;calculate which pixels need to be masked
;to do this change to a new coordinate system (ang)
; 	then work arond the ellipse, each time
;	incrementing angle by some small amount
;	smaller increments for dimmer (smaller) objs.
      FOR ang = 0.0, (2*!PI),increment DO BEGIN
         xedge = sma(d)*cos(theta(d))*cos(ang) - smb(d)*sin(theta(d))*sin(ang)
         yedge = sma(d)*sin(theta(d))*cos(ang) + smb(d)*cos(theta(d))*sin(ang)
;   print, "x, y,xedge, yedge", x,y,xedge,yedge
         
                                ; first get the center pixel
         IF (x(d) GT 0 AND x(d) LT xmax) THEN BEGIN
            IF (y(d) GT 0 AND y(d) LT ymax) THEN BEGIN
               data[x(d),y(d)] = alog10(-1)
                                ;make sure the other pixels are not off the edge
                                ;then mask the little suckers
               IF (xedge+x(d) GT 0 AND xedge+x(d) LT xmax) THEN BEGIN
                  IF (yedge+y(d) GT 0 AND yedge+y(d) LT ymax) THEN BEGIN
                                ; have to do four cases to get everything 
                                ;between x and xegde and y and yedge
                     IF ((xedge+x(d)) LT x(d)) THEN BEGIN
                        If (yedge+y(d) LT y(d)) THEN BEGIN
                           data[xedge+x(d):x(d),yedge+y(d):y(d)] = alog10(-1)
                        ENDIF ELSE BEGIN
                           data[xedge+x(d):x(d),y(d):yedge+y(d)] = alog10(-1)
                        ENDELSE
                     ENDIF ELSE BEGIN
                        If (yedge+y(d) LT y(d)) THEN BEGIN
                           data[x(d):xedge+x(d),yedge+y(d):y(d)] = alog10(-1)
                        ENDIF ELSE BEGIN
                           data[x(d):xedge+x(d),y(d):yedge+y(d)] =alog10(-1)
                        ENDELSE
                     ENDELSE
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         
      ENDFOR                    ; end for each angle
   endfor                       ;end for each detection
   
  ; fits_write, strcompress('/Users/jkrick/Virgo/IRAC/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_myseg_15.fits',/remove_all), data, header

sigmaarr = fltarr(e)
sizearr = fltarr(e)

for n = 1, e do begin 
   size = 1*n
   sumarr = fltarr(n_elements(xcen))
   for i = 0, n_elements(xcen) - 1 do begin 
      goodarea = where(finite(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.]) gt 0, goodcount)
      sumarr(i) = total(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.],/nan) / goodcount
   endfor
   
   sigmaarr(n-1) = stddev(sumarr)
   sizearr(n-1) = (sqrt(goodcount))*0.6
endfor

oplot, sizearr, sigmaarr, thick = 3,  psym = 6

a = where(finite(data) lt 1) 
print, 'twice the area:how many masked pixels',n_elements(a), n_elements(data)

;----------------------------------------------------
;first make the larger masks
;try increasing the mask size of all objects uniformly
   readcol, 'test.cat', NUMBER                 , X_WORLD                , Y_WORLD                , X                , Y                , MAG_AUTO               , FLUX_BEST              , FLUXERR_BEST           , MAG_BEST               , MAGERR_BEST            , BACKGROUND             , FLUX_MAX               , ISOAREA_IMAGE          , ALPHA_J2000            , DELTA_J2000            , A_IMAGE                , B_IMAGE                , THETA_IMAGE            , MU_THRESHOLD           , MU_MAX                 , FLAGS                  , FWHM_IMAGE             , CLASS_STAR             , ELLIPTICITY            

   print, ellipticity(24)

   sma = sqrt(isoarea_image/(!PI*(1-ellipticity)))
   smb = isoarea_image/(!PI*sma)
   sma = 2.0*sma                ; 3.0,2.3 ,1.6                     3.2;3.0*sma;2.6*sma;2*sma
   smb = 2.0*smb                ; 3.0,2.3,1.6                     3.2;3.0*smb;2.6*smb;2*smb
   
;convert theta from degres into radians
   theta = theta_image*(!PI/180)
   increment =  (2*!PI)/500
  
   xmax = 6261                  ;1601;2049;1900;2866;1900;2866;2048;2000;2800;1000
   ymax = 6155                  ;1801;3148;2100;3486;2100;3486;3148;2000;2700;1000
   
   ;for each detection
   for d = long(0), n_elements(number) - 1 do begin

;calculate which pixels need to be masked
;to do this change to a new coordinate system (ang)
; 	then work arond the ellipse, each time
;	incrementing angle by some small amount
;	smaller increments for dimmer (smaller) objs.
      FOR ang = 0.0, (2*!PI),increment DO BEGIN
         xedge = sma(d)*cos(theta(d))*cos(ang) - smb(d)*sin(theta(d))*sin(ang)
         yedge = sma(d)*sin(theta(d))*cos(ang) + smb(d)*cos(theta(d))*sin(ang)
;   print, "x, y,xedge, yedge", x,y,xedge,yedge
         
                                ; first get the center pixel
         IF (x(d) GT 0 AND x(d) LT xmax) THEN BEGIN
            IF (y(d) GT 0 AND y(d) LT ymax) THEN BEGIN
               data[x(d),y(d)] = alog10(-1)
                                ;make sure the other pixels are not off the edge
                                ;then mask the little suckers
               IF (xedge+x(d) GT 0 AND xedge+x(d) LT xmax) THEN BEGIN
                  IF (yedge+y(d) GT 0 AND yedge+y(d) LT ymax) THEN BEGIN
                                ; have to do four cases to get everything 
                                ;between x and xegde and y and yedge
                     IF ((xedge+x(d)) LT x(d)) THEN BEGIN
                        If (yedge+y(d) LT y(d)) THEN BEGIN
                           data[xedge+x(d):x(d),yedge+y(d):y(d)] = alog10(-1)
                        ENDIF ELSE BEGIN
                           data[xedge+x(d):x(d),y(d):yedge+y(d)] = alog10(-1)
                        ENDELSE
                     ENDIF ELSE BEGIN
                        If (yedge+y(d) LT y(d)) THEN BEGIN
                           data[x(d):xedge+x(d),yedge+y(d):y(d)] = alog10(-1)
                        ENDIF ELSE BEGIN
                           data[x(d):xedge+x(d),y(d):yedge+y(d)] =alog10(-1)
                        ENDELSE
                     ENDELSE
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         
      ENDFOR                    ; end for each angle
   endfor                       ;end for each detection
   
;   fits_write, strcompress('/Users/jkrick/Virgo/IRAC/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_myseg_2.fits',/remove_all), data, header

sigmaarr = fltarr(e)
sizearr = fltarr(e)

for n = 1, e do begin 
   size = 1*n
   sumarr = fltarr(n_elements(xcen))
   for i = 0, n_elements(xcen) - 1 do begin 
      goodarea = where(finite(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.]) gt 0, goodcount)
      sumarr(i) = total(data[xcen[i] - size/2.:xcen[i] + size/2., ycen[i] -size/2.:ycen[i] + size/2.],/nan) / goodcount
   endfor
   
   sigmaarr(n-1) = stddev(sumarr)
   sizearr(n-1) = (sqrt(goodcount))*0.6
endfor

oplot, sizearr, sigmaarr, thick = 3,  psym = 5
xyouts, 3.6, 0.00028, 'x', alignment = 0.5
ps_end


print, '2twice the radii'
a = where(finite(data) lt 1) 
print, 'twice the area:how many masked pixels',n_elements(a), n_elements(data)

for c = 0, n_elements(sizearr) - 1 do print, sizearr(c), sigmaarr(c)

end
