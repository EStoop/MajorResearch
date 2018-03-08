pro gyro_darks
!P.multi = [0,1,1]
;read in an AOR of data
;measure centroids.

colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]

ra =     265.12991   
dec=     68.959629

;all the 100s darks in 2011 with this star in at least part of the fov
aor = ['r41608192','r41622272','r41625856','r41679616','r41683200','r41861120','r41864704','r41887488','r41891072','r41929984','r41933568','r41994240','r41997824','r42012160','r42015744','r42039808','r42043392','r42083072','r42086656','r42156800','r42162176','r42553856','r42557696']

;all the 2s darks in 2011 with this star in at least part of the fov
;aor = ['r41607168','r41621248','r41624832','r41678592','r41682176','r41860096','r41863680','r41886464','r41890048','r41928960','r41932544','r41993216','r41996800','r42011136','r42014720','r42038784','r42042368','r42082048','r42085632','r42155776','r42161664','r42552832','r42556672']

for ch = 0, 0 do begin   ;1
   print, 'working on channel', ch + 1
   yarr = fltarr(n_elements(aor)*29)
   xarr = fltarr(n_elements(aor)*29)
   xfwhmarr = fltarr(n_elements(aor)*29)
   yfwhmarr = fltarr(n_elements(aor)*29)
   timearr = fltarr(n_elements(aor)*29)
   c = 0

   pbcdxfwhmarr = fltarr(n_elements(aor))
   pbcdyfwhmarr = fltarr(n_elements(aor))
   pbcdtimearr = fltarr(n_elements(aor))
   p = 0

   for a = 0, n_elements(aor) - 1 do begin
      print, 'working on ', aor(a)
      dir =strcompress( '/Users/jkrick/irac_warm/gyros/darks/t100s/' + string(aor(a)) + '/ch' + string(ch + 1)+'/bcd',/remove_all)
      CD, dir                   ; change directories to the correct AOR directory
                                ;first = 5
      command  =  "find . -name '*cbcd.fits' > ./cbcdlist.txt"
      spawn, command
      command2 =  "find . -name '*cbunc.fits' > ./cbunclist.txt"
      spawn, command2
      
      readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
      readcol,'./cbunclist.txt',uncname, format = 'A', /silent
      
      for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;print, 'working on ', fitsname(i)         
                                ;now where did it really point to?
         fits_read, fitsname(i), im, h
         fits_read, uncname(i), unc, hunc
         get_centroids_for_calstar_jk,im, h, unc, ra, dec,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm, /WARM
        
         ;make a fake fwhm
         adxy, h, ra, dec, x_center, y_center
         ;xfwhm = randomn(seed)
         ;yfwhm = randomn(seed2)
         
                        ;make sure the star is in the FOV
         if x_center gt 8 and y_center gt 8 and x_center lt 248 and y_center lt 248 then begin
            xarr[c] = x_center
            yarr[c] = y_center
            xfwhmarr[c] = xfwhm
            yfwhmarr[c] = yfwhm
            timearr[c] = t  
            c = c + 1
         endif
         
      endfor                    ; for each fits file in the AOR
      
      
   endfor                       ;for each AOR
   
   timearr = timearr[0:c-1]
   xarr = xarr[0:c-1]
   yarr = yarr[0:c-1]
   xfwhmarr = xfwhmarr[0:c-1]
   yfwhmarr = yfwhmarr[0:c-1]
   
;convert timearr into day of year
   timearr = timearr - timearr(0) ; now in seconds since the start
   timearr = timearr / 60. /60./24. ; now in days since the start
   timearr = timearr + 7. ; first darks are on Jan 7, 2011, or day of year 7
   
  before = where(timearr le 134)
  after = where (timearr gt 134)
  print, 'bcd ymean before & after', mean(yfwhmarr[before]), mean(yfwhmarr[after])
  print, 'bcd xmean before & after', mean(xfwhmarr[before]), mean(xfwhmarr[after])

  s = plot(timearr, yfwhmarr,'6r3+',yrange = [1.5, 2.5], xtitle = 'Day of Year 2011',ytitle = 'YFWHM pix')
  line_x = [134, 134]
  line_y = [1.5, 2.5]
  s = polyline(line_x, line_y, thick = 2, color = black,/data,/current)
  ;s = polyline([0,120], [mean(yfwhmarr[before]),mean(yfwhmarr[before])], '--r2', layout = [1,2,1],/data)
  ;s = polyline([145,180], [mean(yfwhmarr[after]),mean(yfwhmarr[after])], '--r2', layout = [1,2,1],/data)

  u = plot(timearr, xfwhmarr,'6r3+',yrange = [1.5, 2.5], xtitle = 'Day of Year 2011',ytitle = 'XFWHM pix')
  u = polyline(line_x, line_y, thick = 2, color = black, /data, /current)
  ;u = polyline([0,120], [mean(xfwhmarr[before]),mean(xfwhmarr[before])], '--r2', layout = [1,2,2])
  ;u = polyline([145,180], [mean(xfwhmarr[after]),mean(xfwhmarr[after])], '--r2', layout = [1,2,2])

    
;-------------------------------------------------------------
;now look at the pbcd's instead of the bcds.
  ;print, 'moving onto pbcd'
   for a = 0, n_elements(aor) - 1 do begin
      ;print, 'pbcd working on ', aor(a)
      dir =strcompress( '/Users/jkrick/irac_warm/gyros/darks/t100s/' + string(aor(a)) + '/ch' + string(ch + 1)+'/pbcd',/remove_all)
      CD, dir                   ; change directories to the correct AOR directory
                              
      command  =  "find . -name '*maic.fits' > ./cbcdlist.txt"
      spawn, command
      command2 =  "find . -name '*munc.fits' > ./cbunclist.txt"
      spawn, command2
       command3 =  "find . -name '*mcov.fits' > ./covlist.txt"
      spawn, command3
      
      readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
      readcol,'./cbunclist.txt',uncname, format = 'A', /silent
      readcol,'./covlist.txt',covname, format = 'A', /silent

      ;print, 'working on ', fitsname
      fits_read, fitsname, im, h
      fits_read, uncname, unc, hunc
      fits_read, covname, covdata, covheader
     get_centroids_for_calstar_jk,im, h, unc, ra, dec,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm, /WARM
        ;make a fake fwhm
     adxy, h, ra, dec, x_center, y_center
         ;xfwhm = randomn(seed)
         ;yfwhm = randomn(seed2)
                                 ;make sure the star is in the FOV
     if covdata(x_center, y_center) gt 1 then begin
        ;print, 'coverage', covdata(x_center, y_center)
        pbcdxfwhmarr[p] = xfwhm
        pbcdyfwhmarr[p] = yfwhm
        pbcdtimearr[p] = t
        p = p + 1
     endif

   endfor

;convert pbcdtimearr into day of year
   pbcdtimearr = pbcdtimearr - pbcdtimearr(0) ; now in seconds since the start
   pbcdtimearr = pbcdtimearr / 60. /60./24. ; now in days since the start
   pbcdtimearr = pbcdtimearr + 7. ; first darks are on Jan 7, 2011, or day of year 7

   pbefore = where(pbcdtimearr le 134)
   pafter = where (pbcdtimearr gt 134)
   print, 'pbcd ymean before & after', mean(pbcdyfwhmarr[pbefore]), mean(pbcdyfwhmarr[pafter])
   print, 'pbcd xmean before & after', mean(pbcdxfwhmarr[pbefore]), mean(pbcdxfwhmarr[pafter])

   s.Select
   s = plot(pbcdtimearr, pbcdyfwhmarr/2.,'6b2o',/overplot,/current, sym_filled = 1)
   s = polyline([0,120], [mean(pbcdyfwhmarr[pbefore])/2.,mean(pbcdyfwhmarr[pbefore])/2.], '--b2',/data,/current)
   s = polyline([145,180], [mean(pbcdyfwhmarr[pafter])/2.,mean(pbcdyfwhmarr[pafter])/2.], '--b2',/data,/current)
   s1 = text(180, 2.1, 'BCD', /data)
   s2 = text(180, 1.6, 'PBCD', /data)
   u.Select
   u = plot(pbcdtimearr, pbcdxfwhmarr/2.,'6b2o',/overplot,  /current, sym_filled = 1)
   u = polyline([0,120], [mean(pbcdxfwhmarr[pbefore])/2.,mean(pbcdxfwhmarr[pbefore])/2.], '--b2',/data)
   u = polyline([145,180], [mean(pbcdxfwhmarr[pafter])/2.,mean(pbcdxfwhmarr[pafter])/2.], '--b2',/data)
  u1 = text(180, 2.1, 'BCD',/data,/current)
   u2 = text(180, 1.6, 'PBCD', /data,/current)

  
 

endfor  ; for each channel

end
