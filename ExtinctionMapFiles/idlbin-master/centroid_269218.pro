pro centroid_269218

  !P.multi = [0,1,1]
;read in an AOR of data
;measure centroids.

  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]
;HD269218
  ra_star =     78.794125 
  dec_star =     -67.275944

  dirname = '/Users/jkrick/irac_warm/r42411776/ch1/bcd'
  ;this is PID 80113
  ;5 - 12 point reuleaux patterns on the same relatively bright star (12s observations).

  CD, dirname                  
  command  =  "find . -name '*cbcd.fits' > ./cbcdlist.txt"
  spawn, command
   
  readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
  yarr = fltarr(n_elements(fitsname) -1)
  xarr = fltarr(n_elements(fitsname)-1)
  timearr = dblarr(n_elements(fitsname)-1)

  for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
     print, 'working on ', fitsname(i)         
     header = headfits(fitsname(i)) ;
     adxy, header, ra_star, dec_star, x, y
     sclk = sxpar(header, 'SCLK_OBS')
 ; print,format='(D26.8)',sclk
     timearr(i-1) = sclk
                                ;now were did it really point to?
     if x gt 8 and x lt 248 and y gt 8 and y lt 248 then begin
        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], ra = ra_star, dec=dec_star 
        xarr(i-1) = x_center
        yarr(i-1) = y_center
        
     endif

   endfor                       ; for each fits file in the AOR
   
  s = plot(xarr, yarr,'6r1+', xtitle = 'X (pixels)',ytitle = 'Y (pixels)', color = 'red', xrange = [0,300], yrange = [0,300])


;now where are the centers supposed to be?
  
  xcenter_reuleaux = [27.0, 37.5, 40.0, 37.5, 27.0, 2.5, -20.0, -40.5, -55.0, -40.5, -20.0, 2.5]       ;in pixels
  ycenter_reuleaux = [-48.0, -24.0, 0.5, 24.5, 48.0, 44.0, 35.5, 19.5, -0.0, -19.0, -35.5, -44.5]       ;in pixels
  
  cluster_x = [0, 76.8, 0, -76.8, 0]       ;in arcsec
  cluster_y = [0, 0, 76.8, 0, -76.8]       ; in arcsec
  
  cluster_x = cluster_x * 1.22
  cluster_y = cluster_y * 1.22
  
  intended_x = fltarr(12*5)
  intended_y = fltarr(12*5)
  
  intended_x = [xcenter_reuleaux+ cluster_x[0], xcenter_reuleaux+ cluster_x[1], xcenter_reuleaux+ cluster_x[2], xcenter_reuleaux+ cluster_x[3], xcenter_reuleaux+ cluster_x[4]]
  intended_y = [ycenter_reuleaux+ cluster_y[0], ycenter_reuleaux+ cluster_y[1], ycenter_reuleaux+ cluster_y[2], ycenter_reuleaux+ cluster_y[3], ycenter_reuleaux+ cluster_y[4]]
  
  intended_x = 128 - intended_x
  intended_y = 128 + intended_y
  
  print, 'intended_y', intended_y
  print, 'observed', yarr
  ;print,format='(D26.8)', timearr ;- timearr(0)
  s = plot( intended_x, intended_y, '6r1+',/overplot, color = 'silver')
  
!P.multi = [0,2,1]
  deltax =  intended_x - xarr
  a = where(xarr ne 0) 
  deltax = deltax[a]
  t = plot(timearr[a] - timearr(0),deltax, '6b1+', xtitle = 'Time (seconds)', ytitle = 'X Delta from intended position')
 
  
  deltay =  intended_y - yarr
  b = where(yarr ne 0) 
  deltay = deltay[b]
  u = plot(timearr[a] - timearr(0),deltay, '6b1+', xtitle = 'Time (seconds)', ytitle = 'Y Delta from intended position')

end
