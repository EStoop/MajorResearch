pro pcrs_latent
  ;ps_open, filename= '/Users/jkrick/irac_warm/pmap/latent.ps',/portrait,/square,/color

  vsym, /polygon, /fill

  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

 
  dir =[ '/Users/jkrick/irac_warm/pmap/0044407296' ]
  for d = 0, 0 do begin  ;look at both aors
     c = 0
     CD, dir[d]                    ; change directories to the correct AOR directory
     command = 'ls IRAC.1*bcd_fp.fits > /Users/jkrick/irac_warm/pmap/bcdlist.txt'
     spawn, command
     
     readcol,'/Users/jkrick/irac_warm/pmap/bcdlist.txt',fitsname, format = 'A', /silent
  
     for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;now were did it really point to?
;        header = headfits(fitsname(i)) ;
;        ra_rqst = sxpar(header, 'RA_RQST')
;        dec_rqst = sxpar(header, 'DEC_RQST')
;        aorlabel = sxpar(header, 'AORLABEL')
        
;        print, fitsname(i), ra_rqst, dec_rqst
;        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_rqst, dec = dec_rqst    
        
;        fluxarr[i] =mean(abcdflux)
;        funcarr[i] = mean(fs)
;        print, i, mean(abcdflux), mean(fs)

        fits_read, fitsname(i), data, header
          ;convert the data from MJy/sr to electrons
        gain = sxpar(header, 'GAIN')
        exptime =  sxpar(header, 'EXPTIME') 
        fluxconv =  sxpar(header, 'FLUXCONV') 
        time = sxpar(header, 'SCLK_OBS')
        sbtoe = gain*exptime/fluxconv
        data = data*sbtoe

        fluxarr =  fltarr(64)
        bkgarr =  fltarr(64)     
        
        for j = 0, 63 do begin
           flux = total(data[13:17, 13:17, j])
           bkg = total(data[5:9, 5:9, j])
           fluxarr[j] = flux
           bkgarr[j] = bkg
        endfor

        print, i, mean(fluxarr), mean(bkgarr)

     endfor                     ; for each fits file in the AOR

  endfor ; for each AOR

  ;ps_close, /noprint,/noid

;part two, look at other AORs taken of this star to figure out if it
;has latents
print, '--------------------------------------'

  dir = '/Users/jkrick/irac_warm/HD165459/r41590784/ch1/bcd' 
     CD, dir                 ; change directories to the correct AOR directory
     command = 'ls *bcd.fits > /Users/jkrick/irac_warm/HD165459/bcdlist.txt'
     spawn, command
     
     readcol,'/Users/jkrick/irac_warm/HD165459/bcdlist.txt',fitsname, format = 'A', /silent
  
     for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        fits_read, fitsname(i), data, header
          ;convert the data from MJy/sr to electrons
        gain = sxpar(header, 'GAIN')
        exptime =  sxpar(header, 'EXPTIME') 
        fluxconv =  sxpar(header, 'FLUXCONV') 
        time = sxpar(header, 'SCLK_OBS')
        sbtoe = gain*exptime/fluxconv
        data = data*sbtoe
        ra_rqst = sxpar(header, 'RA_RQST')
        dec_rqst = sxpar(header, 'DEC_RQST')

        if i lt 9 then begin  ;find the stars center since it is moving in a reuleaux pattern
           get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_rqst, dec = dec_rqst,/silent  
           x = mean(x_center)
           y = mean(y_center)
        endif else begin ; just measure the first star position
           x = 15.
           y = 15.
        endelse
        
        fluxarr =  fltarr(64)
        bkgarr =  fltarr(64)     
        for j = 0, 63 do begin
           flux = total(data[x-2:x+2,y-2:y+2, j])
           bkg = mean(total(data[5:9, 5:9, j])+ total(data[21:25, 21:25, j]))
           fluxarr[j] = flux
           bkgarr[j] = bkg
        endfor
        
        print, i, x, y, mean(fluxarr), mean(bkgarr)
     endfor

 end
  
