pro irac_photvstime
ps_open, filename='/Users/jkrick/irac_warm/cryo_to_warm_mosaic.ps',/portrait,/square,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
colorarr = [redcolor, bluecolor, greencolor, orangecolor, purplecolor, cyancolor, yellowcolor]

ra = [265.00149,265.00754,264.91231,264.95325,264.88594]
dec = [68.995166, 68.991479,69.006478,68.965557,68.975467]

plot, findgen(10), findgen(10), /nodata, xrange = [9E8, 9.4E8], yrange = [.8,1.05], xtitle = 'SCLK_OBS', ytitle = 'Relative flux', ystyle = 1

restore, '/Users/jkrick/iwic/pixel_phase_img_ch1.sav'   

;dirname =[ '/Users/jkrick/irac_warm/r27994368/ch1/pbcd']
;for d = 0, n_elements(dirname) - 1 do begin
;   command = 'ls '+dirname(d) +'/*_bcd.fits > /Users/jkrick/irac_warm/ch1list.txt'
;   spawn, command
;endfor


;read in the list of fits files
; readcol,'/Users/jkrick/irac_warm/ch1list.txt', fitsname, format = 'A', /silent
;if ch eq 1 then readcol,'/Users/jkrick/irac_warm/ch2list.txt', fitsname, format = 'A', /silent

for numstars = 0, n_elements(ra) - 1 do begin

   fitsname = ['/Users/jkrick/irac_warm/r27991552/ch1/pbcd/SPITZER_I1_27991552_0000_3_E5698574_maic.fits', '/Users/jkrick/irac_warm/r27994368/ch1/pbcd/SPITZER_I1_27994368_0000_3_E5701052_maic.fits','/Users/jkrick/irac_warm/r37659136/ch1/pbcd/SPITZER_I1*maic.fits', '/Users/jkrick/irac_warm/r37812480/ch1/pbcd/SPITZER_I1*maic.fits', '/Users/jkrick/irac_warm/r37828096/ch1/pbcd/SPITZER_I1*maic.fits', '/Users/jkrick/irac_warm/r37831936/ch1/pbcd/SPITZER_I1*maic.fits', '/Users/jkrick/irac_warm/r37840896//ch1/pbcd/SPITZER_I1*maic.fits']

;make an object to keep track of numbers

   bdobject = replicate({bdob, xcen:0D,ycen:0D,aperflux:0D, fluxerr:0D, time:0D},n_elements(fitsname))

;for each mosaic
   for f =0, n_elements(fitsname) -1 do begin
      fits_read, fitsname(f), data, header
      print, 'working on ', fitsname(f)
      exptime = sxpar(header, 'EXPTIME')
      flux_conv = sxpar(header, 'FLUXCONV')
      gain = sxpar(header, 'GAIN')
      naxis1 = sxpar(header, 'NAXIS1')
      naxis2 = sxpar(header, 'NAXIS2')
      clock = sxpar(header, 'SCLK_OBS')

      if clock gt 9.3E8 then gain = 3.7
                                   ;convert to electrons
      ;data = data * gain           ;MJy/sr * e/Dn
   
                                ;multiply by exptime
      ;data = data * exptime     ; Mjy/sr *e/Dn *s
   
                                ;divide by fluxconv
     ; data = data / flux_conv   ;= electrons
   
   
      ;convert to microjy
     ; data = data / 34.9837

   
                                ;get rid of NAN's for cntrd and aper!!!
                                ;aper can't handle NAN's inside of the aperture!
      a = where(finite(data) lt 1)
      data[a] = data[a+1]
   
    
      adxy, header, ra[numstars], dec[numstars], x,y
      ;print, 'x, y', x, y
      if x gt 20 and x lt naxis1 - 20 and y gt 20 and y lt naxis2 - 20 then begin
         cntrd, data, x,y,xcen, ycen, 5
         aper, data, xcen, ycen, flux, errap, sky, skyerr, 1, [10], [12,20],/NAN,/flux,/silent
         print, 'aper flux',  flux  ;in electrons

      ;make a correction for array location dependent correction.
      ;need to make a mosaic

      
      ;make a correction for pixel phase 
      ;   xphase = (xcen mod 1) - NINT(xcen mod 1)
      ;  yphase = (ycen mod 1) - NINT(ycen mod 1)
      ;   p = sqrt(xphase^2 + yphase^2)

         ;cryo
      ;   if clock lt 9.3E8 then begin
      ;      correction = 1 + 0.0535*(1/sqrt(2*!Pi) - p)
      ;      corrected_flux = flux / correction
      ;   endif

         ;warm
      ;   if clock gt 9.3E8 then begin
      ;      interpolated_relative_flux = interp2d(relative_flux,x_center,y_center,xphase,yphase,/regular)
      ;      corrected_flux = flux / interpolated_relative_flux
      ;   endif



      ;calculate errors = three terms
         naper = !Pi*10^2        ;number of on target aperture pixels
         nsky = !Pi*20^2 - !Pi*12^2 ;number of sky pixels
         sigmasky = skyerr*naper/sqrt(nsky)

         sigmaskyon =skyerr*sqrt(naper)

         sigmapoisson = sqrt(flux)


         sigmaflux = sqrt(sigmasky^2 + sigmaskyon^2 + sigmapoisson^2)



         ;put flux back into uJy
        ; corrected_flux = corrected_flux *flux_conv
        ; corrected_flux = corrected_flux / exptime
        ; corrected_flux = corrected_flux /gain   ;in Mjy/sr
        ; corrected_flux = corrected_flux / 34.9837   ;in uJy

                                ;need to keep track of the pixel locations and the flux
                                ;can't get too near the edge because the aperture
                                ;photometry need an annulus
         if finite(flux) gt 0 then begin
            bdobject[f].xcen = xcen
            bdobject[f].ycen = ycen
            bdobject[f].aperflux = flux   ;corrected_flux
            bdobject[f].time = clock
            bdobject[f].fluxerr = sigmaflux;*flux_conv/exptime/gain/34.9837
         
         endif else begin
         print, 'not finite', xcen, ycen
                                ;print, 'min', min(data[xcen-20:xcen+20, ycen-20:ycen+20])
;         print, 'ap', mean(data[xcen-10:xcen+10, ycen-10:ycen-10])
;         print, 'sky', mean(data[xcen-20:xcen+20, ycen-20:ycen-20])
         
      endelse
      endif
   
   endfor

   bdobject = bdobject[0:f-1]

   print, bdobject.aperflux, bdobject.fluxerr

;   oploterror, bdobject.time, bdobject.aperflux/ max(bdobject.aperflux), bdobject.fluxerr / max(bdobject.aperflux), color = colors[numstars]
   oplot, bdobject.time, bdobject.aperflux/ max(bdobject.aperflux), color = colorarr[numstars]

endfor


;for the legend
restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'
objnum = findobject(ra, dec)

print, objectnew[objnum].irac1flux
a = fltarr(n_elements(objnum))

legend, psym = [3,3,3,3,3],[string(objectnew[objnum[0]].irac1flux),string(objectnew[objnum[1]].irac1flux),string(objectnew[objnum[2]].irac1flux),string(objectnew[objnum[3]].irac1flux),string(objectnew[objnum[4]].irac1flux)],textcolors = [colorarr[0], colorarr[1], colorarr[2], colorarr[3], colorarr[4]], /top, /right

ps_close, /noprint,/noid

end
