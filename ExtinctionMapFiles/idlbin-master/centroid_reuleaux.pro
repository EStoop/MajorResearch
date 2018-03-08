pro centroid_reuleaux

  !P.multi = [0,1,2]
;read in an AOR of data
;measure centroids.
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
ps_open, filename='/Users/jkrick/irac_warm/reuleaux_36/positions_r42334208.ps',/portrait,/square,/color

  ;colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]
;HD269218
  ra_star =     319.45111;319.48657
  dec_star =    62.674569;62.646679

  dirname = '/Users/jkrick/irac_warm/reuleaux_36/0042334208'

  CD, dirname                  
  command  =  "find . -name 'IRAC.1.*bcd_fp.fits' > ./cbcdlist.txt"
  spawn, command
   
  readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
  yarr = fltarr(144)
  xarr = fltarr(144)
  timearr = dblarr(144)
  for i =1, 144 do begin ;read each cbcd file, find centroid, keep track
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
  
 xcenter_reuleaux = [    8.00,   10.00,   11.33,   12.33,   11.67,   10.00,    8.33,    3.67,   -1.33,   -7.00,  -11.00,  -13.67,  -16.67,  -13.33,  -11.00,   -6.67,   -1.33,    3.67,    9.00,   11.00,   12.33,   12.33,   11.67,    9.00,    5.33,    0.67,   -3.33,   -9.00,  -13.00,  -15.67,  -15.67,  -12.33,   -9.00,   -3.67,    0.67,    5.67]

ycenter_reuleaux = [-15.00,  -9.67,  -5.00,   0.33,   5.00,   9.67,  14.67,  13.33,  12.67,  10.00,   7.33,   4.00,   0.33,  -5.00,  -7.33, -10.33, -12.67, -13.33, -13.00,  -7.67,  -3.00,   2.33,   7.00,  12.67,  14.67,  13.33,  11.67,   9.00,   5.33,   2.00,  -2.67,  -6.00,  -9.33, -11.33, -13.67, -14.33]


;this pattern was observed for 4 repeats at each position

  intended_x = fltarr(n_elements(xcenter_reuleaux)*4)
  intended_y = fltarr(n_elements(xcenter_reuleaux)*4)
  a = 0
  for i = 0, n_elements(intended_x) -1,4 do begin
     intended_x[i] = xcenter_reuleaux[a]
     intended_y[i] = ycenter_reuleaux[a]
    intended_x[i+1] = xcenter_reuleaux[a]
     intended_y[i+1] = ycenter_reuleaux[a]
    intended_x[i+2] = xcenter_reuleaux[a]
     intended_y[i+2] = ycenter_reuleaux[a]
    intended_x[i+3] = xcenter_reuleaux[a]
     intended_y[i+3] = ycenter_reuleaux[a]
     a = a + 1
  endfor

  starx = 82.
  stary = 84.

  intended_x = starx - intended_x
  intended_y = stary + intended_y

  print, 'n_intendted', n_elements(intended_x)
  print, 'n_y', n_elements(xarr)
  print, 'time', n_elements(timearr)
  print, 'intended_x', intended_x
  print, 'observed x', xarr
  ;;print,format='(D26.8)', timearr ;- timearr(0)
  ;s = plot( intended_x, intended_y, '6r1+',/overplot, color = 'silver')
  ; oplot, intended_x, intended_y,  psym = 2, color = bluecolor

;!P.multi = [0,2,1]
  deltax =  intended_x - xarr
  deltax = deltax - deltax[0]
;  a = where(xarr ne 0) 
;  deltax = deltax[a]
 ; t = plot(timearr- timearr(0),deltax, '6b1+', xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position')
  plot, timearr- timearr(0),deltax, psym = 2, xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position', yrange = [-1,1.5] , xrange = [0,2500], charthick = 2, xthick = 2, ythick = 2, thick = 2, title = 'IRU2 Ch1 12s Marengo r42334208'

  
  deltay =  intended_y - yarr
  deltay = deltay - deltay[0]
;  b = where(yarr ne 0) 
;  deltay = deltay[b]
  plot, timearr - timearr(0),deltay, psym =2, xtitle = 'Time (seconds)', ytitle = 'Y Delta from intended position', yrange = [-1,1.5] , xrange = [0,2500], charthick = 2, xthick = 2, ythick = 2, thick = 2

ps_close, /noprint,/noid

end
