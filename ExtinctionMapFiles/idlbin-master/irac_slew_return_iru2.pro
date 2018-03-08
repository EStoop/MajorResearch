pro irac_slew_return_iru2
!P.multi = [0,2,1]

;using the 2 pmap tests done in wk. 509

ps_start, filename='/Users/jkrick/irac_warm/pixel_phase/slew_return_sub_2.ps'
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
blackcolor = FSC_COLOR("Black", !D.Table_Size-3)

ra =     269.74253; 269.730266870728     
dec=     67.794291; 67.7939663925534


sub_aor = ['0044244736','0044244480']


for ch = 1, 1 do begin
   for a = 0,  n_elements(sub_aor) - 1 do begin
      deltax_arr = fltarr(n_elements(sub_aor)*122 )
      deltay_arr = fltarr(n_elements(sub_aor)*122 )
      n = 0
      ni = 0
      print, 'working on ', sub_aor(a)

      ;why are the positions wonky?
      if a eq 0 then begin
         ra =     269.74253     ; 269.730266870728     
         dec=     67.794291     ; 67.7939663925534
      endif
      if a eq 1 then begin
         ra =     269.7274     ; 269.730266870728     
         dec=     67.793567     ; 67.7939663925534
      endif


     ; if ch eq 0 then begin
     ;    dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(sub_aor(a)) 
     ;    CD, dir            ; change directories to the correct AOR directory
      ;   first = 5
      ;   command  =  "ls IRAC.1.*bcd_fp.fits > /Users/jkrick/irac_warm/pmap/bcdlist.txt"
      ;endif
      if ch eq 1 then begin
         dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(sub_aor(a)) 
         CD, dir            ; change directories to the correct AOR directory
         first = 4
         command  =  "ls IRAC.2.*bcd_fp.fits > /Users/jkrick/irac_warm/pmap/bcdlist.txt"
      endif
      
      spawn, command
      
      readcol,'/Users/jkrick/irac_warm/pmap/bcdlist.txt',fitsname, format = 'A', /silent
      yarr = fltarr(n_elements(fitsname))
      xarr = fltarr(n_elements(fitsname))
      x_reqarr = fltarr(n_elements(fitsname))
      y_reqarr = fltarr(n_elements(fitsname))
      c = 0
;since they bounce back and forth, need to know which frame actually has the star on it.

      for i =first,n_elements(fitsname) - 2, 2 do begin ;read each cbcd file, find centroid, keep track
         ;print, 'working on ', fitsname(i)
         if i eq first then header_master = headfits(fitsname(i)) ;
         header = headfits(fitsname(i))
                                ;figure out where the reqeusted pointing was
         ra_req= sxpar(header, 'RA_RQST')
         dec_req = sxpar(header, 'DEC_RQST')
         adxy, header_master, ra_req, dec_req, x_req, y_req
         x_reqarr[c] = x_req
         y_reqarr[c] = y_req
         
                                ;now were did it really point to?
         get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], RA=ra, DEC=dec ,/silent
         xarr[c] = mean(x_center)
         yarr[c] = mean(y_center)
         c = c + 1
      endfor ; for each fits file in the AOR
      
      xarr = xarr[0:c-1]
      yarr = yarr[0:c-1]
      x_reqarr = x_reqarr[0:c-1]
      y_reqarr = y_reqarr[0:c-1]
      
       ; from one position to the next, what was the intended shift, and the measured shift? 
      for ni = 1, n_elements(xarr) - 1 do begin
         deltax_arr[n] = abs(xarr[ni] - xarr[ni-1]) - abs( x_reqarr[ni] - x_reqarr[ni - 1])
         print, 'x measured, intended', abs(xarr[ni] - xarr[ni-1]), abs(x_reqarr[ni] - x_reqarr[ni - 1])
         deltay_arr[n] = abs(yarr[ni] - yarr[ni-1]) - abs(y_reqarr[ni] - y_reqarr[ni - 1])
         print, 'y measured, intended', abs(yarr[ni] - yarr[ni-1]), abs(y_reqarr[ni] - y_reqarr[ni - 1])
        n = n + 1
      endfor

   deltax_arr = deltax_arr[0:n-1]
   deltay_arr = deltay_arr[0:n-1]


   plothist, deltax_arr, xxhist, xyhist, bin = 0.005, xrange = [-0.2, 0.2], xtitle = 'Delta x', ytitle = 'Number', title = 'Ch' + string(ch + 1),/noplot
;fit with a gaussian
   start = [0.01,0.02, 500]
   noise = fltarr(n_elements(xxhist)) + 1.0
   result= MPFITFUN('mygauss',xxhist,xyhist, noise, start)    
   yaxis = result(2)/sqrt(2.*!Pi) * exp(-0.5*((xxhist - result(0))/result(1))^2.)
   yaxis = yaxis / max(yaxis) ; normalize
;   plot,xxhist, yaxis , xrange = [-0.1, 0.1], xtitle = 'Actual Shift - Intended Shift (pixels)',  xthick = 3, ythick = 3, ytitle = 'Number', title = string(sub_aor(a))
   plot,xxhist, xyhist , xrange = [-0.1, 0.1], xtitle = 'Actual Shift - Intended Shift (pixels)',  xthick = 3, ythick = 3, ytitle = 'Number', title = string(sub_aor(a)), psym = 10, yrange = [0,40]
   ;xyouts, 0.1, 0.8, result(1)
   
   
   plothist, deltay_arr, yxhist, yyhist, bin = 0.01, xrange = [-0.2, 0.2], xtitle = 'Delta y', ytitle = 'Number', title = 'Ch' + string(ch + 1),/noplot
;fit with a gaussian
  start = [0.01,0.02, 500]
  noise = fltarr(n_elements(yxhist)) + 1.0
  resulty= MPFITFUN('mygauss',yxhist,yyhist, noise, start)    
  yaxis =  resulty(2)/sqrt(2.*!Pi) * exp(-0.5*((yxhist - resulty(0))/resulty(1))^2.)
  yaxis = yaxis / max(yaxis)    ; normalize
;  oplot, yxhist,yaxis,linestyle = 2
   oplot, yxhist,yyhist,linestyle = 2, psym = 10
 ;xyouts, 0.1,1, result(1)

al_legend, ['x', 'y'], linestyle = [0,2],/top,/left

   endfor ;for each AOR
   
   
endfor                          ; for each channel

;-----------------------
;and for the legend


ps_end

end
