pro irac_slew_return_2
!P.multi = [0,1,1]
;read in an AOR of data
;measure centroids.
;this is NpM1p67.0536
;****these need to be corrected
ps_start, filename='/Users/jkrick/irac_warm/pixel_phase/slew_return_sub_2.ps'
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
blackcolor = FSC_COLOR("Black", !D.Table_Size-3)

ra =     269.730266870728     
dec=     67.7939663925534

full_aor = ['38933504','38933760','38934016','38934272','38934528','38934784','38935040','38935296','38935552','38935808','39113984','39115264','39114240'];,'39114496','39114752','39115008','39115520','39115776','39116032','39116288']
;full_aor = ['38933504','38933760']

sub_aor = ['38936064','38936320','38936576','38936832','39112960','39113472','39113216','39113728']

for ch = 1, 1 do begin
   deltax_arr = fltarr(n_elements(sub_aor)*122 )
   deltay_arr = fltarr(n_elements(sub_aor)*122 )
   n = 0
   for a = 0, n_elements(sub_aor) - 1 do begin
      ni = 0
      print, 'working on ', sub_aor(a)
      if ch eq 0 then begin
         dir = '/Users/jkrick/irac_warm/pixel_phase/r' + string(sub_aor(a)) + '/ch1/bcd'
         CD, dir            ; change directories to the correct AOR directory
         first = 4
      endif
      if ch eq 1 then begin
         dir = '/Users/jkrick/irac_warm/pixel_phase/r' + string(sub_aor(a)) + '/ch2/bcd'
         CD, dir            ; change directories to the correct AOR directory
         first = 5
      endif
      
      command  =  "find . -name '*bcd.fits' > ./cbcdlist.txt"
      spawn, command
      
      readcol,'./cbcdlist.txt',fitsname, format = 'A', /silent
      yarr = fltarr(n_elements(fitsname))
      xarr = fltarr(n_elements(fitsname))
      x_reqarr = fltarr(n_elements(fitsname))
      y_reqarr = fltarr(n_elements(fitsname))
      c = 0
;for channel 1, I know that frame 1 has the star, and every odd frame.

      for i =first,n_elements(fitsname) - 1, 2 do begin ;read each cbcd file, find centroid, keep track
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
         xarr[c] = x_center
         yarr[c] = y_center
         
                                ;some diagnostics
                                ;print, i, ' ',x_center - x_req, ' ',x_center, ' ', x_req, ' ', y_req, format = '(I10, A, F10.6, A, F10.6, A, F10.6,A, F10.6)'
         
         c = c + 1
      endfor ; for each fits file in the AOR
      
      xarr = xarr[0:c-1]
      yarr = yarr[0:c-1]
      x_reqarr = x_reqarr[0:c-1]
      y_reqarr = y_reqarr[0:c-1]
      
      
      
      for ni = 1, n_elements(xarr) - 1 do begin
         deltax_arr[n] = abs(xarr[ni] - xarr[ni-1]) - abs( x_reqarr[ni] - x_reqarr[ni - 1])
         print, 'measured, intended', abs(xarr[ni] - xarr[ni-1]), abs(x_reqarr[ni] - x_reqarr[ni - 1])
         deltay_arr[n] = abs(yarr[ni] - yarr[ni-1]) - abs(y_reqarr[ni] - y_reqarr[ni - 1])
         n = n + 1
      endfor


   endfor ;for each AOR
   
   deltax_arr = deltax_arr[0:n-1]
   deltay_arr = deltay_arr[0:n-1]
   plothist, deltax_arr, xxhist, xyhist, bin = 0.01, xrange = [-0.2, 0.2], xtitle = 'Delta x', ytitle = 'Number', title = 'Ch' + string(ch + 1),/noplot
;fit with a gaussian
   start = [0.01,0.02, 500]
   noise = fltarr(n_elements(xxhist)) + 1.0
   result= MPFITFUN('mygauss',xxhist,xyhist, noise, start)    
   yaxis = result(2)/sqrt(2.*!Pi) * exp(-0.5*((xxhist - result(0))/result(1))^2.)
   yaxis = yaxis / max(yaxis) ; normalize
   plot,xxhist, yaxis , xrange = [-1.0, 1.0], xtitle = 'Actual Position - Intended Position (pixels)',  xthick = 3, ythick = 3, ytitle = 'Number'
   ;xyouts, 0.1, 0.8, result(1)
   
   
   plothist, deltay_arr, yxhist, yyhist, bin = 0.01, xrange = [-0.2, 0.2], xtitle = 'Delta y', ytitle = 'Number', title = 'Ch' + string(ch + 1),/noplot
;fit with a gaussian
  start = [0.01,0.02, 500]
  noise = fltarr(n_elements(yxhist)) + 1.0
  result= MPFITFUN('mygauss',yxhist,yyhist, noise, start)    
  yaxis =  result(2)/sqrt(2.*!Pi) * exp(-0.5*((yxhist - result(0))/result(1))^2.)
  yaxis = yaxis / max(yaxis)    ; normalize
  oplot, yxhist-0.02,yaxis,linestyle = 2
  ;xyouts, 0.1,1, result(1)
   
endfor                          ; for each channel

;------------------------------
;now add the random slew information

pixel_dist = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
xhist = [3.0,5.0,3.0,10.0,9.0,15.0,28.0,42.0,55.0,85.0,136.0,144.0,143.0,135.0,108.0,68.0,65.0,42.0,28.0,16.0,14.0]
yhist = [12.0,17.0,44.0,40.0,53.0,64.0,85.0,81.0,129.0,136.0,136.0,121.0,78.0,64.0,37.0,18.0,16.0,11.0,6.0,4.0,2.0]

;plot, pixel_dist, yhist, xtitle = 'Distance from Desired Position in Pixels', ytitle = 'Number', /nodata, thick = 3, charthick = 3, xthick = 3, ythick = 3
;oplot, pixel_dist+.108, yhist, color = redcolor, thick = 3
;oplot,  pixel_dist - .18, xhist, color = bluecolor, thick = 3

start = [0.1,0.2, 5000.]
noise = fltarr(n_elements(yhist))
noise[*] = 1                                                    ;equally weight the values

;x
resultx= MPFITFUN('mygauss',pixel_dist,xhist, noise, start)    ;./quiet    

y = (resultx(2))/sqrt(2.*!Pi) * exp(-0.5*((pixel_dist - (resultx(0)))/(resultx(1)))^2.)
y = y / max(y)
oplot, pixel_dist - 0.18, y,  thick  = 3, color  = redcolor

;y
resulty= MPFITFUN('mygauss',pixel_dist,yhist, noise, start)    ;./quiet    
y = (resulty(2))/sqrt(2.*!Pi) * exp(-0.5*((pixel_dist - (resulty(0)))/(resulty(1)))^2.)
y = y / max(y)
oplot, pixel_dist +.108, y , thick  = 3, color = redcolor,linestyle =  2 




;-----------------------
;and for the legend
;legend, ['0.31', '0.35','0.04', '0.02'], linestyle = [2,1,2,1], color = [redcolor, redcolor, blackcolor, blackcolor],/top,/right

ps_end

end
