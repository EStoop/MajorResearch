pro centroid_reuleaux_pregyro

  !P.multi = [0,1,2]
;read in an AOR of data
;measure centroids.
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
;ps_open, filename='/Users/jkrick/irac_warm/reuleaux_36/positions_r41976064.ps',/portrait,/square,/color
ps_open, filename='/Users/jkrick/irac_warm/reuleaux_36/positions_r39440640.ps',/portrait,/square,/color

  ;colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]
;HD269218
  ra_star = 138.2329; 53.104028
  dec_star = -1.0201506;  -27.803929

  ;dirname = '/Users/jkrick/irac_warm/reuleaux_36/r41976064/ch1/bcd/'
  dirname = '/Users/jkrick/irac_warm/reuleaux_36/r39440640/ch1/bcd/'

  CD, dirname                  
  command  =  "find . -name '*cbcd.fits' > ./cbcdlist.txt"
  spawn, command
   
  readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
  yarr = fltarr(72)
  xarr = fltarr(72)
  timearr = dblarr(72)
  for i =1, 72 do begin ;read each cbcd file, find centroid, keep track
     ;print, 'working on ', fitsname(i)         
     header = headfits(fitsname(i)) ;
     adxy, header, ra_star, dec_star, x, y
     sclk = sxpar(header, 'SCLK_OBS')
 ; print,format='(D26.8)',sclk
     timearr(i-1) = sclk
                                ;now were did it really point to?
     if x gt 8 and x lt 248 and y gt 8 and y lt 248 then begin
        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], ra = ra_star, dec=dec_star ,/silent
        xarr(i-1) = x_center
        yarr(i-1) = y_center
        
     endif

   endfor                       ; for each fits file in the AOR

  ;s = plot(xarr, yarr,'6r1+', xtitle = 'X (pixels)',ytitle = 'Y (pixels)', color = 'red', xrange = [140,180], yrange = [180,220])
 ; plot, xarr, yarr, psym = 2, xtitle = 'X (pixels)',ytitle = 'Y (pixels)';, color = redcolor;, xrange = [140,180], yrange = [180,220]

;now where are the centers supposed to be?
  ;medium pattern

; xcenter_reuleaux = [  16.00,  18.00,  20.33,  22.33,  23.67,  24.00,  24.33,  24.67,  23.67,  22.00,  20.00,  18.33,  16.33,  11.67,   6.00,   1.33,  -3.33,  -7.33, -13.00, -17.00, -20.67, -24.67, -27.33, -31.00, -33.67, -30.33, -27.33, -25.00, -21.00, -16.67, -12.67,  -7.33,  -4.00,   1.33,   6.67,  11.67]


; ycenter_reuleaux = [-29.00,-24.67,-20.00,-14.67,-10.00, -5.33,  0.67,  5.33, 10.67, 14.00, 19.33, 24.00, 28.33, 28.00, 27.67, 26.67, 25.33, 23.67, 21.00, 18.33, 15.00, 11.33,  8.00,  4.67,  0.67, -4.67, -8.33,-12.00,-15.67,-19.00,-21.67,-24.00,-25.33,-26.33,-27.67,-28.33]

 xcenter_reuleaux = [   33.00 ,   37.00 ,   41.33 ,   44.33 ,   46.67 ,   48.00 ,   48.33 ,   48.67 ,   46.67 ,   44.00 ,   41.00 ,   37.33 ,   33.33 ,   23.67 ,   13.00 ,    3.33 ,  -6.330 ,  -15.33 ,  -25.00 ,  -33.00 ,  -40.67 ,  -48.67 ,  -54.33 ,  -62.00 ,  -66.67 ,  -61.33 ,  -54.33 ,  -49.00 ,  -41.00 ,  -32.67 ,  -24.67 ,  -15.33 ,  -7.000 ,    3.33 ,   13.67 ,   23.67 ]


 ycenter_reuleaux = [-58.00,-48.67,-40.00,-29.67,-20.00,-10.33,  0.67, 10.33, 20.67, 29.00, 39.33, 48.00, 57.33, 57.00, 55.67, 53.67, 50.33, 46.67, 42.00, 36.33, 30.00, 23.33, 16.00,  8.67,  0.67, -8.67,-16.33,-24.00,-30.67,-37.00,-42.67,-47.00,-50.33,-53.33,-55.67,-57.33]



;this pattern was observed for 2 repeats at each position
;ignore the second observation
; a = 0
; good_xarr = fltarr(n_elements(xarr) / 2)
; good_yarr = fltarr(n_elements(yarr) / 2)
; for i = 2, n_elements(xarr) - 1,2 do begin
;    print, 'xarr[i],i, a', xarr[i],i, a
;    good_xarr[a] = xarr[i]
;    good_yarr[a] = yarr[i]
;    a = a + 1
;endfor
 good_xarr = xarr
 good_yarr = yarr

  starx = 45.;170
  stary = 86.;140

  intended_x = starx - xcenter_reuleaux
  intended_y = stary + ycenter_reuleaux

  print, 'n_intendted', n_elements(intended_x)
  print, 'n_y', n_elements(xarr)
  print, 'time', n_elements(timearr)
  print, 'intended_x', intended_x
  print, 'observed x', xarr
  ;;print,format='(D26.8)', timearr ;- timearr(0)
  ;s = plot( intended_x, intended_y, '6r1+',/overplot, color = 'silver')
  ; oplot, intended_x, intended_y,  psym = 2, color = bluecolor

;!P.multi = [0,2,1]
  deltax =  intended_x - good_xarr
  deltax = deltax - deltax[0]
;  a = where(xarr ne 0) 
;  deltax = deltax[a]
 ; t = plot(timearr- timearr(0),deltax, '6b1+', xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position')
  plot, timearr- timearr(0),deltax, psym = 2, xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position', yrange = [-1,1.5] , xrange = [0,2500], charthick = 2, xthick = 2, ythick = 2, thick = 2, title = 'IRU1 Ch1 30s Cooray r39440640'

  
  deltay =  intended_y - good_yarr
  deltay = deltay - deltay[0]
;  b = where(yarr ne 0) 
;  deltay = deltay[b]
  plot, timearr - timearr(0),deltay, psym =2, xtitle = 'Time (seconds)', ytitle = 'Y Delta from intended position', yrange = [-1,1.5], xrange = [0,2500], charthick = 2, xthick = 2, ythick = 2, thick = 2


ps_close, /noprint,/noid

end
