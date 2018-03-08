pro bright_star_test
!P.multi=[0,1,2]

ps_open, filename='/Users/jkrick/iwic/brightstartest_0.ps',/portrait,/square,/color
;ps_start, filename='/Users/jkrick/iwic/brightstartest.ps'
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

;measure the latents on the very bright star tests

;mag 1 test is in IWICDAA
;dir_name = '/Users/jkrick/iwic/iwicDAA/IRAC020700'
;aor_list = 'aor_vbright_mag1'
;mag 0 test is in IWICDAI
;dir_name = '/Users/jkrick/iwic/iwicDAI/IRAC021300'
;aor_list = 'aor_vbright_mag0'
dir_name = '/Users/jkrick/iwic/iwicDAS/IRAC021900'
aor_list = 'aor_vbright_magsuper'

;organize all the filenames 
filenames, dir_name, aor_list


flux_conv = [.0858219,.127683]  ;mag1
;flux_conv = [.0858219,.128090] ; mag 0 (uses different biases in ch2)

exptime = [2,12,100]

;load the darks in units of DN
fits_read, '/Users/jkrick/iwic/iwicDAI/IRAC021300/darks/ch1.skydrk100s.fits', dark_ch1_100s, darkhead_ch1_100s
fits_read, '/Users/jkrick/iwic/iwicDAI/IRAC021300/darks/ch1.skydrk12s.fits', dark_ch1_12s, darkhead_ch1_12s
fits_read, '/Users/jkrick/iwic/iwicDAI/IRAC021300/darks/ch1.skydrk2s.fits', dark_ch1_2s, darkhead_ch1_2s
fits_read, '/Users/jkrick/iwic/iwicDAI/IRAC021300/darks/ch2.skydrk100s.fits', dark_ch2_100s, darkhead_ch2_100s
fits_read, '/Users/jkrick/iwic/iwicDAI/IRAC021300/darks/ch2.skydrk12s.fits', dark_ch2_12s, darkhead_ch2_12s
fits_read, '/Users/jkrick/iwic/iwicDAI/IRAC021300/darks/ch2.skydrk2s.fits', dark_ch2_2s, darkhead_ch2_2s

;plotimage, xrange=[1,255],yrange=[1,255], bytscl(dark_ch1_100s[*,*,0]) ,/preserve_aspect, /noaxes, ncolors=8
;need to flip the darks about the y-axis

dark_ch1_100s = reverse( dark_ch1_100s[*,*,0], 2)
dark_ch1_12s = reverse( dark_ch1_12s[*,*,0], 2)
dark_ch1_2s = reverse( dark_ch1_2s[*,*,0], 2)
dark_ch2_100s = reverse( dark_ch2_100s[*,*,0], 2)
dark_ch2_12s = reverse( dark_ch2_12s[*,*,0], 2)
dark_ch2_2s = reverse( dark_ch2_2s[*,*,0], 2)

;what is the x, y location of the bright star?
;figure this out more precisely later

xcen1 = 11
ycen1 =  38
xcen2 = 247
ycen2 =  37


bcdimage = 'test'
readcol, dir_name +'/'+ aor_list, aorname, format="A",/silent

for ch = 0,   1 do begin        ; 1 do begin       ;for ch1 and ch2
   count = 0
   obstime = dblarr(1000)
   flux = fltarr(1000)
   indpix = fltarr(1000)
   indpix1 = fltarr(1000)
   indpix2 = fltarr(1000)
   indpix3 = fltarr(1000)
   indpix4 = fltarr(1000)
   indpix5 = fltarr(1000)
   if ch eq 0 then begin
      xcen = xcen1
      ycen = ycen1
   endif else begin
      xcen = xcen2
      ycen = ycen2
   endelse

   for pos = 0,   n_elements(aorname)  -1  do begin ; 
      
      for exp =  0, n_elements(exptime) - 1 do begin
         print, '-----------'
         print, 'working on ch, exp ', ch+ 1, ' ', exptime(exp)
         
         ;choose the appropriate dark, and put it in DN/s
         if ch eq 0 and exp eq 0 then dark = dark_ch1_2s / float(exptime(exp))
         if ch eq 0 and exp eq 1 then dark = dark_ch1_12s / float(exptime(exp))
         if ch eq 0 and exp eq 2 then dark = dark_ch1_100s/ float(exptime(exp))
         if ch eq 1 and exp eq 0 then dark = dark_ch2_2s/ float(exptime(exp))
         if ch eq 1 and exp eq 1 then dark = dark_ch2_12s/ float(exptime(exp))
         if ch eq 1 and exp eq 2 then dark = dark_ch2_100s/ float(exptime(exp))


         if ch eq 0 then begin
            bcdlist =  strcompress(dir_name + '/bcd/'+ aorname(pos)+ '/bcdlist_ch1_' + string(exptime(exp))+'s.txt',/remove_all)
         endif else begin
            bcdlist =  strcompress(dir_name + '/bcd/'+ aorname(pos)+ '/bcdlist_ch2_' + string(exptime(exp))+'s.txt',/remove_all)
         endelse
         print, 'bcdlist', bcdlist

         ;it is ok to have empty files
         openr, emptylun, bcdlist, error = err, /get_lun
         while  ~ EOF(emptylun) do begin

            readf, emptylun, bcdimage

            fits_read, dir_name+ '/bcd/'+ aorname(pos)+'/'+bcdimage, bcddata, bcdheader ;'/Users/jkrick/nutella/IRAC/iwic210/'
            indpix(count) = bcddata[xcen, ycen]
            indpix1(count) = bcddata[xcen+2, ycen]
            indpix2(count) = bcddata[xcen+4, ycen]
            indpix3(count) = bcddata[xcen+6, ycen]
            indpix4(count) = bcddata[202,144]
            if ch eq 0 then  indpix5(count) = bcddata[212,148]  ;go after an unaffected pixel
            if ch eq 1 then indpix5(count) = bcddata[78,232]
                                ;correct the central pixel to a
                                ;ballpark idea of what it should be
                                ;without saturation
            if bcdimage eq './IRAC.1.0036894208.0003.0000.1.bcd_fp.fits' then bcddata[xcen, ycen] = 1.5E6
            if bcdimage eq './IRAC.2.0036894464.0003.0000.1.bcd_fp.fits' then bcddata[xcen, ycen] = 1.3E6

            ;convert to DN/s
            bcddata = bcddata / flux_conv(ch)
            
            ;subtract away the dark
            bcddata = bcddata - dark
            
            ;convert back into MJy/sr
            bcddata = bcddata * flux_conv(ch)

            newname = strmid(bcdimage, 2, 31) +  'bcd_dark.fits'
            newname = strcompress(dir_name + '/bcd/' +  aorname(pos)+'/'+newname)
            print, 'writing ', newname
            fits_write, newname, bcddata, bcdheader

            ;measure aperture photometry and keep track of the UTCS_OBS
            obstime(count) =    fxpar (bcdheader, 'UTCS_OBS')
            print, 'ch, exp,obstime', ch, exptime(exp), obstime(count)

            object = bcddata[xcen -4:xcen +4,ycen -4:ycen +4]
            nobject = n_elements(object)
            if ch eq 0 then begin
               sky = total(bcddata[xcen :xcen +70, ycen :ycen +70])- total(bcddata[xcen :xcen +60,ycen :ycen +60])  
               nannulus = n_elements(bcddata[xcen :xcen +70, ycen :ycen +70]) - n_elements(bcddata[xcen :xcen +60,ycen :ycen +60])    
            endif else begin
               sky = total(bcddata[xcen -70:xcen , ycen :ycen +70])- total(bcddata[xcen - 60 :xcen ,ycen :ycen +60])  
               nannulus = n_elements(bcddata[xcen -70 :xcen , ycen :ycen +70]) - n_elements(bcddata[xcen - 60 :xcen ,ycen :ycen +60])    
            endelse

            sky = sky / nannulus ;mean sky flux
            flux(count) =  total( object)  - sky*nobject

            count = count + 1
         endwhile
         close,emptylun
         free_lun, emptylun
     
      endfor                    ; for each exptime

   endfor                       ;for each aor
   obstime = obstime[0:count - 1]
   flux = flux[0:count - 1]
   print, flux
   plot, obstime - min(obstime), flux/max(flux), psym  = 2, xtitle = 'time in seconds ', $
         ytitle = 'fractional flux', yrange=[0.000001,1], /ylog, title = 'Ch '+ string(ch + 1) + '   Mag 0 RR UMi', thick = 3, $
         xthick=3, ythick = 3, charthick = 3;, xrange=[1,1000]

 ; plot, findgen(20), findgen(20), /nodata
   vsym, /polygon, /fill


   indpix = indpix[0:count - 1]
   indpix1 = indpix1[0:count - 1]
   indpix2 = indpix2[0:count - 1]
   indpix3 = indpix3[0:count - 1]
   indpix4 = indpix4[0:count - 1]
   indpix5 = indpix5[0:count - 1]


;   plot, obstime - min(obstime),indpix- indpix5, psym = 8, symsize = 0.8, xtitle = 'time in seconds ',$
;         ytitle = ' Delta flux',  title = 'Ch '+ string(ch + 1) + '   Mag 0 RR UMi', thick = 3, $
;         xthick=3, ythick = 3, charthick = 3, yrange=[-4,4];, xrange=[1,1000]
 
;   oplot, obstime - min(obstime),indpix1 - indpix5, psym = 8, symsize = 0.8, color = redcolor
;   oplot, obstime - min(obstime),indpix2- indpix5, psym = 8, symsize = 0.8, color =orangecolor
;   oplot, obstime - min(obstime),indpix3- indpix5, psym = 8, symsize = 0.8, color = greencolor
   oplot, obstime - min(obstime),indpix4- indpix5, psym = 8, symsize = 0.8, color = bluecolor
;   oplot, obstime - min(obstime),indpix5;, psym = 8, symsize = 0.8, color = bluecolor
   print, 'min(obstime)', min(obstime)
endfor                          ;for each channel

ps_close, /noprint,/noid
;ps_end, /png

end
