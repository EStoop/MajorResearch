pro pcrs_inc_c
!P.multi = [0,1,1]
;read in an AOR of data
;measure centroids.

colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY'    ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU', ' PINK',  'PLUM',  ' POWDER_BLUE' ]

ra =     269.72751   
dec=     67.79463

;second test in april
sub_aor = ['0042781696','0042781184', '0042781440',  '0042780160',  '0042779904',  '0042780416', '0042780672',  '0042783488', '0042781952',  '0042780928']


deltaxarr =  fltarr(n_elements(sub_aor))
deltayarr =  fltarr(n_elements(sub_aor))
for a = 0, n_elements(sub_aor) - 1 do begin
   ni = 0
   print, 'working on ', sub_aor(a)
   dir = '/Users/jkrick/irac_warm/pcrs_inc_c/' + string(sub_aor(a)) 
   

   CD, dir                      ; change directories to the correct AOR directory
                                ;first = 5
   ;command  =  "find . -name '*bcd_fp.fits' > ./cbcdlist.txt"
   command  =  "find . -name '*bcd_fp.fits' > ./cbcdlist.txt"
   spawn, command
   
   readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
   yarr = fltarr(n_elements(fitsname)*64)
   xarr = fltarr(n_elements(fitsname)*64)
   c = 0

   for i =0, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
      get_centroids_xy,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], XCEN = 15, YCEN = 15 ,/silent
      header = headfits(fitsname(i)) ;
      if i eq 0 then begin
         xarr = x_center
         yarr = y_center
                                ;plot, xarr, yarr,xrange = [14.0,15.5], yrange = [14.5,16]
                                ;print, 'time', time
      endif else begin
         xarr = [xarr, x_center]
         yarr = [yarr, y_center]
      endelse
      
      deltaxarr[a] = mean(x_center) ;- 15.0
      deltayarr[a] = mean(y_center) ;- 15.0
      ;print,'mean(x_center, ycenter', mean(x_center), mean(y_center)
     ; if i eq 1 then oplot, x_center, y_center, color = bluecolor
     ;       if i eq 2 then oplot, x_center, y_center, color = greencolor
     ; if i eq 3 then begin
      ;   oplot, x_center, y_center, color = yellowcolor
      ;   print, 'xcenter', x_center
       ;  print, 'ycenter', y_center
      ;endif

      ;if i eq 4 then oplot, x_center, y_center, color = cyancolor

   endfor                       ; for each fits file in the AOR
   
   timearr = findgen(n_elements(xarr)) * 0.4

   if a eq 0 then begin
      ;t = plot(xarr, yarr, '6r1*',xrange = [14.0,15.5], yrange = [14.0,16.0], xtitle = 'X pix',ytitle = 'Y pix', color = colorarr[a])
      s = plot(timearr, yarr,'6r1+',yrange = [14.0, 16.0], xrange = [0, 450], xtitle = 'Time(seconds)',ytitle = 'Y pix', color = colorarr[a])
   endif

   if a gt 0 then begin
      ;t = plot( xarr, yarr, '6r1*',/overplot, color = colorarr[a])
      s = plot( timearr, yarr, '6r1+',/overplot, color = colorarr[a])
   endif

endfor                          ;for each AOR

;for a = 0, n_elements(sub_aor) -1 do begin
;   ttext = text( deltaxarr[a], deltayarr[a],sub_aor[a],/data, color = colorarr[a], font_size = 10)
;endfor
;!p.multi = [0,1,2]
;plothist, deltaxarr, xhist, yhist, bin=0.04, xtitle = 'xcenter'
;start = [14.78,0.1, 6]
;noise = fltarr(n_elements(xhist))
;noise = noise + 1.
;result= MPFITFUN('mygauss',xhist,yhist, noise, start)    ;ICL
;oplot, xhist,  result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), thick = 3,linestyle = 2;, color = colors.green

;plothist, deltayarr, xhist, yhist, bin=0.04, xtitle = 'ycenter'

;d = plot(deltaxarr, t0arr/60./60., 'b2*',xtitle = 'delta x from center (pix)', ytitle = 'time since s2pcal (hrs)', yrange = [0,12], xrange = [-2.0,0.5])
;d2 = plot(deltayarr, t0arr/60./60., 'b2*',xtitle = 'delta y from center (pix)', ytitle = 'time since s2pcal (hrs)', yrange = [0,12], xrange = [-0.5, 1.0])
end
