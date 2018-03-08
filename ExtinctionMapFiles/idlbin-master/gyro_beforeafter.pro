pro gyro_beforeafter
  !P.multi = [0,1,1]
  ps_open, filename='/Users/jkrick/irac_warm/gyros/agnreverb/ycenter_time.ps',/portrait,/square,/color

                                ;keep track of how well the telescope
                                ;is executing a dither pattern on agnreverb data
                                ;before and after gyro failure started happening

;from http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/dither/311_cycling_med.tab
  dithernum = findgen(311)
  cyclingx =[  32.0,   4.5,  11.0, -34.5,  35.0,   8.5, -50.0,  10.5,  10.0, -52.5,  35.0, -35.5,  74.0, -46.5,  16.0,  30.5, -35.0,  34.5,  18.0,  20.5,   3.0,  13.5, -29.0,  -4.5,   0.0, -62.5,  23.0,  13.5, -15.0, -49.5,  -0.0, -31.5,  -3.0,  39.5,   8.0,   7.5,  22.0,  -1.5,   4.0, -17.5, -47.0, -17.5,  20.0,  32.5,   8.0,  -6.5, -37.0,  13.5,  33.0,  -9.5,  53.0,-100.5,   1.0, -23.5,   8.0,  -7.5,  30.0, -24.5, -22.0,  46.5,  25.0,   0.5,   1.0,  79.5, -18.0,  -0.5,  39.0,  36.5,   7.0, -37.5, -15.0, -38.5,  21.0,  51.5,  24.0, -21.5,  -5.0,  49.5,  -7.0, -43.5,  11.0,  66.5, -21.0,  31.5, -14.0,  -2.5,   1.0,  53.5, -24.0,  -3.5,  33.0,  31.5, -42.0,  10.5,  44.0,  -5.5, -10.0, -26.5, -18.0, -29.5,  20.0,   1.5,  26.0,  -4.5, -35.0,  76.5,  37.0,  19.5, -19.0,  -1.5,   0.0,   1.5,  28.0, -34.5,   8.0, -15.5,  11.0,  -9.5,   7.0,  47.5, -30.0,  29.5, -24.0,  20.5,  25.0,   5.5, -46.0,   5.5,  24.0,  30.5, -25.0,  43.5, -26.0,  21.5,  37.0,  30.5, -14.0,   6.5,   5.0,  13.5, -16.0,  25.5, -30.0, -71.5, -18.0,  -5.5,  -5.0,  26.5,   7.0,  47.5, -41.0, -24.5,  11.0,  -9.5,  21.0, -22.5,  12.0, -36.5, -24.0, -12.5, -10.0, -48.5,  -9.0,  14.5,  -2.0,  -3.5,  29.0,  61.5,  35.0,  -9.5,  21.0,  26.5,  48.0,  -8.5,  33.0,  -5.5,  23.0, -58.5, -27.0,  22.5, -35.0, -27.5,  15.0,   7.5, -50.0,  -9.5, -24.0,  -2.5,  -0.0, -45.5, -10.0,  -5.5,  38.0, -62.5,  42.0,  38.5,  54.0,  19.5,  -6.0,  74.5,  12.0, -23.5,  53.0, -31.5,  56.0, -12.5, -42.0,  -3.5,  35.0,  21.5,  33.0, -24.5,  23.0, -10.5,  29.0,  11.5,   0.0,  11.5,  35.0,  57.5,  19.0,  62.5, -21.0,  20.5,  17.0,  25.5, -24.0, -17.5, -34.0,  -0.5,  -8.0,  37.5,  16.0, -35.5,   6.0, -25.5, -23.0,  23.5,  72.0,  28.5,  -3.0, -41.5,  55.0,  84.5,  23.0, -21.5,   4.0, -11.5,  10.0,  32.5,  70.0, -10.5,  -9.0,  23.5, -19.0,  -6.5, -81.0,  71.5,  46.0, -64.5, -28.0,  -6.5,  35.0, -22.5, -14.0,  54.5, -77.0,  27.5,  -4.0, -20.5,  14.0,  31.5,  22.0,  67.5,  -9.0,  57.5,  39.0,  37.5,  -4.0, -23.5,  -0.0,   4.5, -33.0, -12.5, -27.0,   8.5,   3.0, -17.5,  36.0, -56.5,  23.0, -31.5, -28.0,  12.5,  -7.0,  50.5,  50.0, -16.5, -12.0,  -7.5,  38.0, -36.5,  83.0,  35.5,  -3.0, -38.5,  15.0,   3.5,  -8.0,   4.5, -21.0]
  cyclingy = [ -40.0,   2.0,  41.5,  87.5, -29.0,  29.0,   4.5,  19.5, -32.0,  10.0, -12.5, -63.5, -17.0,  62.0,  11.5, -59.5,  21.0,   0.0,  18.5, -14.5,  25.0,   0.0,   5.5, -15.5,   6.0,   0.0,  57.5, -16.5,  33.0,  -4.0,  -2.5, -15.5, -10.0,  -8.0,   3.5,  22.5,  41.0,  69.0, -30.5,  23.5, -32.0, -27.0, -11.5, -10.5, -46.0,  22.0, -16.5,  34.5,  63.0, -86.0,  14.5, -12.5,  46.0, -48.0, -38.5, -11.5,  -1.0, -34.0,  25.5,  -9.5, -18.0, -33.0, -30.5,   8.5,  42.0,  18.0,  -3.5,   8.5,  38.0,  22.0,  -0.5,   6.5,  16.0,  13.0,   6.5,  13.5,  45.0,  -2.0, -14.5,   0.5,  13.0,  12.0, -32.5, -11.5,  31.0,  62.0,  28.5,-106.5,   5.0,   4.0,   5.5,   0.5, -48.0,  -3.0,  32.5, -44.5,  -5.0, -11.0, -11.5,   8.5,  32.0,  -5.0,  21.5,  48.5,   0.0,  -1.0,  48.5, -52.5, -33.0,  30.0,  25.5,  19.5,  42.0, -98.0,  53.5,  36.5,  -2.0,  60.0,  36.5,  39.5,  -9.0,  27.0, -62.5,  36.5,  -5.0, -29.0, -22.5,  14.5,  13.0, -23.0,  16.5, -14.5,  10.0,  23.0,  16.5, -38.5,  -6.0,   4.0,  20.5,  -5.5,  39.0, -26.0, -32.5,  42.5,  64.0,  -3.0,  29.5,  -4.5, -13.0,  77.0, -23.5, -17.5,  19.0, -33.0, -42.5,  15.5, -22.0,  51.0,  -1.5,  -7.5,  13.0,  -4.0,  21.5,  29.5,  71.0,  20.0,  50.5, -26.5, -44.0,  28.0,  -5.5,  25.5,  -6.0, -42.0,  24.5,  -1.5, -61.0,  -1.0,  41.5, -14.5, -15.0,  -1.0,   7.5,  51.5,  -7.0,  19.0,  50.5, -92.5,  16.0, -26.0,  -9.5, -77.5, -25.0,  34.0,  55.5, -18.5, -14.0,  15.0, -28.5, -28.5, -14.0,  38.0,  24.5, -17.5,  35.0, -18.0, -11.5, -19.5,  31.0, -33.0, -18.5, -41.5,  61.0, -11.0,  -7.5,  38.5,  -6.0,  50.0,  14.5, -18.5,  -9.0,   8.0, -34.5,  14.5,   4.0,  11.0,  14.5,  22.5,  12.0,  13.0, -22.5,  34.5, -32.0,  43.0,  13.5,  31.5,  26.0, -50.0,   4.5,  -5.5,  31.0, -55.0, -15.5, -40.5, -14.0, -15.0,  -1.5,   3.5,  54.0,  46.0,  33.5,  51.5, -45.0, -34.0,  -5.5,  -6.5, -29.0, -48.0, -47.5,  -4.5,  32.0, -60.0, -17.5,  17.5,  -3.0, -19.0,  20.5,  34.5, -30.0, -62.0, -24.5,  12.5, -27.0,  -3.0, -25.5,   6.5,  26.0,   5.0,   5.5,   3.5,  36.0, -25.0,   2.5,  35.5, -18.0, -38.0, -20.5,  36.5, -39.0,  57.0,  47.5,  30.5, -33.0,  31.0,  42.5,  -3.5, -20.0, -12.0,  -3.5, -61.5,  19.0,  28.0,   9.5, -20.5, -13.0, -14.0,  20.5, -37.5, -18.0, -22.0,  38.5]

  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]
;HD269218
  ra_star =    219.06933; 219.09231 
  dec_star =     58.798641; 58.794373
  
  aor = ['r40185344','r40185344','r40185600','r40185856','r40186112','r40186368','r40186624','r40186880','r40187136','r40187392','r40187648','r40187904','r40188160','r40188416','r40188672','r40188928','r40189184','r40189440','r40189696','r40189952','r40190208']

  for r = 0, n_elements(aor) -1 do begin
     dirname = strcompress('/Users/jkrick/irac_warm/gyros/agnreverb/'+aor(r)+'/ch1/bcd',/remove_all)
     CD, dirname                  
     command  =  "find . -name '*cbcd.fits' > ./cbcdlist.txt"
     spawn, command
     command2  =  "find . -name '*cbunc.fits' > ./cbunclist.txt"
     spawn, command2
     
     readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'./cbunclist.txt',uncname, format = 'A', /silent
     yarr = fltarr(n_elements(fitsname)/2)
     xarr = fltarr(n_elements(fitsname)/2)
     timearr = dblarr(n_elements(fitsname)/2)
     c = 0
     for i =1, n_elements(fitsname),2 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
        header = headfits(fitsname(i)) ;
        adxy, header, ra_star, dec_star, x, y
        sclk = sxpar(header, 'SCLK_OBS')
                                ; print,format='(D26.8)',sclk
        timearr(c) = sclk
                                ;now were did it really point to?
        if x gt 8 and x lt 248 and y gt 8 and y lt 248 then begin
           get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], ra = ra_star, dec=dec_star ,/silent
         fits_read, fitsname(i), im, h
         fits_read, uncname(i), unc, hunc
         ;get_centroids_for_calstar_jk,im, h, unc, ra_star, dec_star,  t, dt, hjd, xft, x3, y3, $
         ;                          x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
          ;                         x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
         ;                          xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
         ;                          xfwhm, yfwhm, /WARM
 
           xarr(c) =  x_center
           yarr(c) =  y_center
          ; print, 'x,y fwhm', xfwhm, yfwhm
           c = c + 1
        endif
        
     endfor                     ; for each fits file in the AOR
     
                                ;s = plot(xarr, yarr,'6r1+', xtitle = 'X (pixels)',ytitle = 'Y (pixels)', color = 'red', xrange = [0,255], yrange = [0,255])
    ;if r eq 0 then begin
    ;   plot, xarr, yarr,psym = 1, xtitle = 'X (pixels)',ytitle = 'Y (pixels)', xrange = [0,250], yrange = [0,250], xcharsize = 2, ycharsize = 2
    ;endif else begin
    ;   oplot, xarr, yarr, psym =1
    ;endelse

;now where are the centers supposed to be?
     
     intended_x = cyclingx[233:252]
     intended_y = cyclingy[233:252]
     
;is the 'star' off center
     intended_x = 141 - intended_x ;130
     intended_y = 198 + intended_y ;127
     
     ;print, 'n_intendted', n_elements(intended_x)
     ;print, 'n_y', n_elements(xarr)
     ;print, 'time', n_elements(timearr)
     ;print, 'intended_y', intended_y
     ;print, 'observed y', yarr
     ;;print,format='(D26.8)', timearr ;- timearr(0)
                                ;s = plot( intended_x, intended_y, '6r1+',/overplot, color = 'silver')
     ;oplot, intended_x, intended_y, psym = 2
     
     deltax =  intended_x - xarr
;  a = where(xarr ne 0) 
;  deltax = deltax[a]
                                ;t = plot(timearr- timearr(0),deltax, '6b1+', xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position')
     ;if r eq 0 then begin
     ;   plot, timearr- timearr(0),deltax, xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position', xcharsize = 2, ycharsize = 2, yrange = [-1, 1]
     ;endif
     ;if r le 10 and r gt 0 then begin
        ;oplot, timearr - timearr(0),deltax
        ;print, 'should be making solid lines'
    ; endif
    ; if r gt 10 then begin
     ;   oplot, timearr - timearr(0),deltax, linestyle = 1
     ;endif

     deltay =  intended_y - yarr
     ;print, 'delta y', deltay - deltay(0)
     ;print, 'time', timearr - timearr(0)
     ;print ,'rrrrrrrrrrrrrr', r
     ;print, 'aor', aor(r)
;  b = where(yarr ne 0) 
;  deltay = deltay[b]
                                ;u = plot(timearr - timearr(0),deltay, '6b1+', xtitle = 'Time (seconds)', ytitle = 'Y Delta from intended position')
     if r eq 0 then begin 
        plot, timearr - timearr(0),deltay - deltay(0),  xtitle = 'Time (seconds)', ytitle = 'Y Delta from intended position', xcharsize = 2, ycharsize = 2, yrange = [-1, 1], charthick = 3, thick = 3, xthick = 3, ythick = 2
     endif 
     if r le 10 and r gt 0 then begin
        oplot, timearr - timearr(0),deltay- deltay(0), thick = 2

     endif
     if r gt 10 then begin
        oplot, timearr - timearr(0),deltay- deltay(0), linestyle = 1, thick = 3
     endif


  endfor  ;for each AOR
ps_close, /noprint,/noid

  
end
