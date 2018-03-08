pro pmap

 
  dir =[ '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/0044244736'  ,'/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/0044244480']
  for d = 0, 0 do begin  ;look at both aors
     c = 0
     CD, dir[d]                    ; change directories to the correct AOR directory
     command = 'ls IRAC.2*bcd_fp.fits > /Users/jkrick/irac_warm/pmap/bcdlist.txt'
     spawn, command
     
     readcol,'/Users/jkrick/irac_warm/pmap/bcdlist.txt',fitsname, format = 'A', /silent
     yarr =   fltarr(n_elements(fitsname)*64)
     xarr =  fltarr(n_elements(fitsname)*64)
     timearr =  fltarr(n_elements(fitsname)*64)
     c = 0
     
    if d eq 0 then begin
         ra =     269.74253     ; 269.730266870728     
         dec=     67.794291     ; 67.7939663925534
      endif
      if d eq 1 then begin
         ra =     269.7274     ; 269.730266870728     
         dec=     67.793567     ; 67.7939663925534
      endif
      

     for i =4, n_elements(fitsname) - 1,2 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra, dec = dec, /silent
        
        header = headfits(fitsname(i)) ;
        sclk = sxpar(header, 'SCLK_OBS')
        
        for j = 0, n_elements(x_center) - 1 do begin
           xarr[c] = x_center[j]
           yarr[c] = y_center[j]
           timearr[c] = sclk + ( 0.1*j )
           ;if y_center[n_elements(y_center) - 1] - y_center[1] gt 0.3 then print, 'inside',  y_center[j]
           c = c + 1
        endfor
        
     endfor                     ; for each fits file in the AOR
     
     timearr = timearr[0:c-1]
     xarr = xarr[0:c-1]
     yarr = yarr[0:c-1]

     timearr = timearr - timearr[0] 
  
     if d eq 0 then begin
        s = plot(timearr, yarr,'6r1+' ,xtitle = 'Time(s)',ytitle = 'Y pix', yrange = [13.0, 16.0], xrange = [4000,6000])
        t = plot(timearr, xarr,'6r1+' ,xtitle = 'Time(s)',ytitle = 'X pix', yrange = [14.0, 17.0], xrange = [4000,6000])
     endif
     if d eq 1 then begin
        s = plot(timearr, yarr,'6r1+' ,xtitle = 'Time(s)',ytitle = 'Y pix' , yrange = [10.0,22.0], xrange = [0,7000])
        t = plot(timearr, xarr,'6r1+' ,xtitle = 'Time(s)',ytitle = 'X pix', yrange = [10.0, 22.0], xrange = [0,7000])
     endif
     

;figure out where the jumps are.
for count = 1, n_elements(xarr) - 3 do begin

   ax = where(xarr[count] - xarr[count +1] gt 0.4 , axcount)
   if axcount gt 0 then print, 'found an x break', xarr[count], xarr[count + 1], timearr[count]

   ay = where(yarr[count+1] - yarr[count ] gt 0.4 , aycount)
   if aycount gt 0 then print, 'found a y break', yarr[count], yarr[count + 1], timearr[count]

endfor


  endfor ; for each AOR


 end
  
