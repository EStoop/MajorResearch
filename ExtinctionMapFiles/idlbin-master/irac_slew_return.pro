pro irac_slew_return
!P.multi = [0,2,2]
;read in an AOR of data
;measure centroids.
;this is NpM1p67.0536
;****these need to be corrected
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
ps_start, filename='/Users/jkrick/irac_warm/pixel_phase/slew_return_sub.ps'

ra =     269.730266870728     
dec=     67.7939663925534

full_aor = ['38933504','38933760','38934016','38934272','38934528','38934784','38935040','38935296','38935552','38935808','39113984','39115264','39114240'];,'39114496','39114752','39115008','39115520','39115776','39116032','39116288']
;full_aor = ['38933504','38933760']

sub_aor = ['38936064','38936320','38936576','38936832','39112960','39113472','39113216','39113728']

for ch = 0, 1 do begin
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
   plothist, deltax_arr, xxhist, xyhist, bin = 0.01, xrange = [-0.2, 0.2], xtitle = 'Delta x', ytitle = 'Number', title = 'Ch' + string(ch + 1)
;fit with a gaussian
   start = [0.01,0.02, 500]
   noise = fltarr(n_elements(xxhist)) + 1.0
   result= MPFITFUN('mygauss',xxhist,xyhist, noise, start)    
   oplot,xxhist,  result(2)/sqrt(2.*!Pi) * exp(-0.5*((xxhist - result(0))/result(1))^2.), linestyle = 2
   xyouts, 0.1, 100, result(1)
   
   
   plothist, deltay_arr, yxhist, yyhist, bin = 0.01, xrange = [-0.2, 0.2], xtitle = 'Delta y', ytitle = 'Number', title = 'Ch' + string(ch + 1)
;fit with a gaussian
  start = [0.01,0.02, 500]
  noise = fltarr(n_elements(yxhist)) + 1.0
  result= MPFITFUN('mygauss',yxhist,yyhist, noise, start)    
  oplot, yxhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((yxhist - result(0))/result(1))^2.),linestyle = 2
  xyouts, 0.1,200, result(1)
   
endfor                          ; for each channel

ps_end

end
