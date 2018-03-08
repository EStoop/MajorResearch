pro pmap_test
  ps_open, filename= '/Users/jkrick/irac_warm/pmap/five_diamonds.ps',/portrait,/square,/color

  vsym, /polygon, /fill

  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

  ra =     270.627943             
  dec=     58.627237         

  dir =[ '/Users/jkrick/irac_warm/pmap/0044407040' ]
  for d = 0, 0 do begin  ;look at both aors
     c = 0
     CD, dir[d]                    ; change directories to the correct AOR directory
     command = 'ls IRAC.1*bcd_fp.fits > /Users/jkrick/irac_warm/pmap/bcdlist.txt'
     spawn, command
     
     readcol,'/Users/jkrick/irac_warm/pmap/bcdlist.txt',fitsname, format = 'A', /silent
     yarr =   fltarr(n_elements(fitsname)*64)
     xarr =  fltarr(n_elements(fitsname)*64)
     timearr =  fltarr(n_elements(fitsname)*64)
     x_meanarr =  fltarr(n_elements(fitsname)*64)
     y_meanarr =  fltarr(n_elements(fitsname)*64)
 c = 0
     
 
     for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;now were did it really point to?
;      if i eq 1 then begin
;         header = headfits(fitsname(i)) ;
;         ra_ref = sxpar(header, 'RA_REF')
;         dec_ref = sxpar(header, 'DEC_REF')
;         aorlabel = sxpar(header, 'AORLABEL')
;      endif
      
      
;      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent     
;      print, i, mean(x_center), mean(y_center)
;     if i eq 1 then begin 
;        plot, x_center, y_center, psym = 8, xrange =[14.5, 15.5], yrange = [14.5, 15.5], ystyle =1, xstyle = 1, title = aorlabel, symsize = 0.5
;     endif else begin
;        oplot, x_center, y_center, psym = 8, symsize = 0.5
;     endelse
     
;     x_meanarr[i] = mean(x_center)
;     y_meanarr[i] = mean(y_center)



     endfor                     ; for each fits file in the AOR
       
     xyouts,x_meanarr, y_meanarr, 'x', alignment = 0.5, color = redcolor


  endfor ; for each AOR

  ps_close, /noprint,/noid


 end
  
