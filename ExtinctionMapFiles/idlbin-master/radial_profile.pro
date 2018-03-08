;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;February 2002
;Jessica Krick
;
;This program creates a radial profile of a star in a fits image.  
; first it blocks out the bad column and diffraction spikes, if given,
; then it displays the image on the screen.  The radial distance from the
; center, and the counts at that distance are stored in an array. counts 
; above and below a threshold are not included in the array.  The array
; is binned into 1 pixel radii, and the counts are averaged within the bin.
; The counts are then normalized and plotted.
; 
;
;input: filename, xcenter, ycenter, location of bad columns, and saturated stars
;output: fielname.prof, idlprof.ps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro radial_profile

close, /all		;close all files = tidiness
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
!P.multi = [0,1,2]
;ps_open, filename='/Users/jkrick/Virgo/IRAC/s18p14/psf.ps',/portrait,/square,/color
ps_start, filename= '/Users/jkrick/Virgo/IRAC/psf_masked.ps',/nomatch

!P.thick = 3
!P.charthick = 3

for ch = 0, 1 do begin
  
;declare variables
   if ch eq 0 then begin 
      xcenter = ulong(1890.5)   ; x coord of center of the star
      ycenter = ulong(5325.6)   ; y coord of center of the star
      filename = "/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/mosaic_bkgd"
      print, 'ch1**************'
   endif
   
   if ch eq 1 then begin
      xcenter = ulong(1222.8)   ; x coord of center of the star
      ycenter = ulong(5338.5)   ; y coord of center of the star
      filename = "/Users/jkrick/Virgo/IRAC/ch2/ch2_Combine-mosaic/mosaic_bkgd"
      print, 'ch2**************'

   endif

   x = ulong(0)			; x pixel coord
   y = ulong(0)			; y pixel coord
   r = float(0.0)               ; radius from center in pixels
   p = float(0.0)               ; count at a particular radius
   max = 400.0
   imagefile = filename + '.fits'
   datafile = filename + '.prof'
; read in image
   FITS_READ, imagefile, data, header

;try increasing the mask size of all objects uniformly
   CD, strcompress('/Users/jkrick/Virgo/irac/s18p14/ch' + string(ch+ 1) + '/ch' + string(ch+1) + '_Combine-mosaic',/remove_all)

   readcol, 'withoutstar.cat', NUMBER                 , X_WORLD                , Y_WORLD                , X                , Y                , MAG_AUTO               , FLUX_BEST              , FLUXERR_BEST           , MAG_BEST               , MAGERR_BEST            , BACKGROUND             , FLUX_MAX               , ISOAREA_IMAGE          , ALPHA_J2000            , DELTA_J2000            , A_IMAGE                , B_IMAGE                , THETA_IMAGE            , MU_THRESHOLD           , MU_MAX                 , FLAGS                  , FWHM_IMAGE             , CLASS_STAR             , ELLIPTICITY            

   print, ellipticity(24)

   sma = sqrt(isoarea_image/(!PI*(1-ellipticity)))
   smb = isoarea_image/(!PI*sma)
   sma = 2.*sma                ; 3.0,2.3 ,1.6                     3.2;3.0*sma;2.6*sma;2*sma
   smb = 2.*smb                ; 3.0,2.3,1.6                     3.2;3.0*smb;2.6*smb;2*smb
  
;convert theta from degres into radians
   theta = theta_image*(!PI/180)
   increment =  (2*!PI)/500
 
   xmax = 6261                  ;1601;2049;1900;2866;1900;2866;2048;2000;2800;1000
   ymax = 5975                  ;1801;3148;2100;3486;2100;3486;3148;2000;2700;1000
  
   ;for each detection
   for d = long(0),  n_elements(number) - 1 do begin
;      if y(d) gt 5330 and y(d) lt 5347 and x(d) gt 1217 and x(d) lt 1230 then print, x(d), y(d), flux_best(d)

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
  
   fits_write, strcompress('/Users/jkrick/Virgo/IRAC/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_starseg.fits',/remove_all), data, header
   fits_read, strcompress('/Users/jkrick/Virgo/IRAC/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_starseg.fits',/remove_all), data, header

;---------------------------------------------------------------

;get a real centroid
   ;xyad, header, xcenter, ycenter, ra1, dec1
   ;print, 'ra1 ',  ra1, dec1
;get_centroids,imagefile, t, dt, xcenter, ycenter, flux, xs, ys, fs, b, /WARM, /APER, APRAD = 3, RA=ra1, DEC=dec1,/silent
;   cntrd, data, xcenter, ycenter, x, y, 3.0
;   print, 'x,y', x, y
;   xcenter = x
;   ycenter = y
; make an array of radii from the central pixel, and counts for the 
; corresponding radii

;before looping, check to make sure that the edges of the box are
;not off the image, whoops!
; !!!!!!!!!!!!!!need to add one for greater than 2048
   IF (xcenter - max LT 0) THEN BEGIN
      max = float(xcenter - 0)
      print, "radius range is not allowed in x, reseting max = ", max
   ENDIF
   
   IF (ycenter - max LT 0) THEN BEGIN
      max = float(ycenter - 0)
      print, "radius range is not allowed in y, reseting max = ", max
   ENDIF
   
   radius = FLTARR((4*max)^2)	; array of radius (float)
   counts = FLTARR((4*max)^2)	; array of counts (float)
   
; two for loops get all x and y in a square = 2*max around the center
;it then outputs the radius and counts to a file and into two arrays
   
   j = ulong64(0)               ;counter
   increment = 1                ;sampling distance
   FOR x = xcenter - max, xcenter + max, increment DO BEGIN
      FOR y = ycenter -max, ycenter + max, increment DO BEGIN
                                ;print, x, y
         s = float(xcenter - x)
         t = float(ycenter - y)
         r = float(sqrt(s^(2)+(t^2))) ;pythagoras
         p = data[x,y]
         IF (p GT -1 AND P LT 33000) THEN BEGIN ;block the bad data pts
                                ;printf,lun,r,p			;output to datafile
            counts(j) = p
            radius(j) = r
            j= j + 1
         ENDIF
                                ;print, x, y    
      ENDFOR 
   ENDFOR
   print, "j = ", j
   
;shorten the arrays to be only as long as need be
   radius = radius(0:j-1)
   counts = counts(0:j-1)
   
;sort the radius array, and then apply that sorting order to the counts
;this way both arrays get sorted appropriately
   sortindex = Sort(radius)
   sortedradius = radius[sortindex]
   sortedcounts = counts[sortindex]
   
;try binning with histogram instead of manually.
   
   h=histogram(sortedradius,BINSIZE=1, REVERSE_INDICES=ri)
   print, 'n_element(h)', n_elements(h)
   meanarr=fltarr(n_elements(h))
   radarr=fltarr(n_elements(h))
   for j=0L,n_elements(h) -1 do begin
      if ri[j+1] gt ri[j] then begin
         meanarr[j]=median(sortedcounts[ri[ri[j]:ri[j+1]-1]])
         radarr[j]=mean(sortedradius[ri[ri[j]:ri[j+1]-1]])
      endif
      
   endfor
   
   
;plot
   if ch eq 0 then begin
      plot, radarr*0.6, meanarr/meanarr(1), xtitle = 'Radius in Arcsec', ytitle = 'Normalized Flux', /xlog, /ylog, xrange = [1, 1000], yrange = [1E-8,1], thick = 3, charthick =3, xthick = 3, ythick = 3
      p1 = !P & x1 = !X & y1 = !Y
      legend, ['Ch1', 'Ch2'], linestyle = [0,2],/bottom
      a = findgen(1000) 
      oplot, a, ((a)^(-2.5)) /20, linestyle = 1, thick = 3

   endif else begin
      !P = p1 & !X = x1 & !Y = y1
      oplot,  radarr*0.6, meanarr/meanarr(1), linestyle = 2
   endelse
   ;a = findgen(1000) + 3
;oplot, a, ((a)^(-3)) * 10^(3.8), color = bluecolor
   
   
;try encircled energy
;integrate under the psf profile
   ee = fltarr(n_elements(meanarr))
   ee2 = fltarr(n_elements(meanarr))
   for r = 1, n_elements(radarr) - 1 do begin
      
      ee(r) = int_tabulated(radarr(0:r)*0.6, meanarr(0:r),/double)
      ee2(r) = tsum(radarr(0:r)*0.6, meanarr(0:r))
      
   endfor
   
   
   if ch eq 0 then begin
      plot, radarr*0.6, ee2/ee2(n_elements(ee2)-1),/xlog,  xrange = [1, 1000], yrange = [0.99, 1.0], xtitle =  'Radius in Arcsec', ytitle = 'Encircled Energy' , thick = 3, charthick =3, xthick = 3, ythick = 3
      p2 = !P & x2 = !X & y2 = !Y
   endif else begin
      !P = p2 & !X = x2 & !Y = y2
      oplot, radarr*0.6, ee2/ee2(n_elements(ee2)-1), linestyle = 2
   endelse



endfor
;ps_close, /noprint,/noid
ps_end
END

