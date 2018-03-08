pro run_latent, dir_name, aor_bright
 ;run_latent, '/Users/jkrick/iwic/iwicDA/IRAC019000', 'aor_t37r35b60_50_bright'
!p.multi = [0, 2, 2]

ps_open, filename='/Users/jkrick/iwic/iwicDA/IRAC019000/latent.ps',/portrait,/square,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
graycolor = FSC_COLOR("gray", !D.Table_Size-9)
browncolor = FSC_COLOR("brown", !D.Table_Size-10)

colorname = [ graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor,graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor,graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor ]

;track the locations on the chip on which the bright sources are, and
;see if their aperture fluxes vary over time.

;read in the correct data

;organize input files

readcol, dir_name + '/'+ aor_bright, brightaor, format="A",/silent

;two star positions
;ra = [271.015833, 270.93958333]
;dec = [66.928333, 66.93416667 ]

ra = [271.296500,271.354292]
dec = [66.945500,66.913333]
gain = [3.3,3.7]

for ca = 0 , n_elements(brightaor) - 2 do begin

   command1 = ' find '+dir_name+'/bcd/'+brightaor(ca)+' -name "*.fits" > ' + dir_name + '/bcd/'+brightaor(ca)+'/files_test.list'
   command2 = 'grep IRAC.1 < '+dir_name+'/bcd/'+brightaor(ca)+'/files_test.list  | grep -v IRAC.2 | grep bcd_fp.fits > '+dir_name+'/bcd/'+brightaor(ca)+'/ch1_bcd_test.list'
   command3 = 'grep IRAC.2 < '+dir_name+'/bcd/'+brightaor(ca)+'/files_test.list  | grep bcd_fp.fits > '+dir_name+'/bcd/'+brightaor(ca)+'/ch2_bcd_test.list'

   a = [command1, command2, command3]
   for i = 0, n_elements(a) -1 do spawn, a(i)

   readcol, dir_name+'/bcd/'+brightaor(ca)+'/ch1_bcd_test.list', bcdname_ch1, format="A",/silent
   readcol, dir_name+'/bcd/'+brightaor(ca)+'/ch2_bcd_test.list', bcdname_ch2, format="A",/silent
   
  
;-----------------------------------------------------   
;forchannel 1
   k = 0
   bcdname = bcdname_ch1

   flux_star1_pos1 = fltarr(20)
   flux_star1_pos2 = fltarr(12)
   flux_star2_pos1 = fltarr(20)
   flux_star2_pos2 = fltarr(12)
   
   
   imagenum = 4
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos1(0) = flux(0)
   flux_star2_pos1(0) = flux(1)
   
   getrot, header, rot, cdelt              ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                     ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0)    ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])           ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)               ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1)    ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])           ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)               ;Within 5 arc seconds
   
   for imagenum = 5, 22 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;aper is being a shit and won't work for the positions without
;a star on them

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos1(imagenum -4) = flux1
      flux_star2_pos1(imagenum -4) = flux2
   endfor
   ;------------------

   imagenum = 12
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos2(0) = flux(0)
   flux_star2_pos2(0) = flux(1)
   
   getrot, header, rot, cdelt                  ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                         ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0) ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])             ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)                 ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1) ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])             ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)                 ;Within 5 arc seconds

   for imagenum = 13, 22 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos2(imagenum -12) = flux1
      flux_star2_pos2(imagenum -12) = flux2
   endfor


plot, findgen(20), flux_star1_pos1 / max(flux_star1_pos1), /ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch1 star1 0.4s'
oplot, findgen(12), flux_star1_pos2 / max(flux_star1_pos2), thick=3, color=bluecolor

plot, findgen(20), flux_star2_pos1/ max(flux_star2_pos1),/ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch1 star2 0.4s'
oplot, findgen(12), flux_star2_pos2/ max(flux_star2_pos2), thick=3, color=bluecolor
;-----------------------------------------------------   
;forchannel 2
   k = 1
   bcdname = bcdname_ch2

   flux_star1_pos1 = fltarr(24)
   flux_star1_pos2 = fltarr(16)
   flux_star1_pos3 = fltarr(8)
   flux_star2_pos1 = fltarr(24)
   flux_star2_pos2 = fltarr(16)
   flux_star2_pos3 = fltarr(8)
   
   
   imagenum = 0
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos1(0) = flux(0)
   flux_star2_pos1(0) = flux(1)
   
   getrot, header, rot, cdelt              ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                     ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0)    ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])           ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)               ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1)    ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])           ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)               ;Within 5 arc seconds
   
   for imagenum = 1, 22 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
      
;aper is being a shit and won't work for the positions without
;a star on them

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos1(imagenum ) = flux1
      flux_star2_pos1(imagenum ) = flux2
   endfor
   ;------------------

   imagenum = 8
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos2(0) = flux(0)
   flux_star2_pos2(0) = flux(1)
   
   getrot, header, rot, cdelt                  ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                         ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0) ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])             ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)                 ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1) ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])             ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)                 ;Within 5 arc seconds

   for imagenum = 9, 22 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos2(imagenum -8) = flux1
      flux_star2_pos2(imagenum -8) = flux2
   endfor

   ;------------------

   imagenum = 16
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos2(0) = flux(0)
   flux_star2_pos2(0) = flux(1)
   
   getrot, header, rot, cdelt                  ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                         ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0) ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])             ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)                 ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1) ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])             ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)                 ;Within 5 arc seconds

   for imagenum = 17, 22 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos3(imagenum -16) = flux1
      flux_star2_pos3(imagenum -16) = flux2
   endfor

plot, findgen(24), flux_star1_pos1 / max(flux_star1_pos1), /ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch2 star1 0.4s'
oplot, findgen(16), flux_star1_pos2 / max(flux_star1_pos2), thick=3, color=bluecolor
oplot, findgen(8), flux_star1_pos3 / max(flux_star1_pos3), thick=3, color=redcolor

plot, findgen(24), flux_star2_pos1/ max(flux_star2_pos1),/ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch2 star2 0.4s'
oplot, findgen(16), flux_star2_pos2/ max(flux_star2_pos2), thick=3, color=bluecolor
oplot, findgen(8), flux_star2_pos3/ max(flux_star2_pos3), thick=3, color=redcolor

endfor

;-----------------------------------------------------   
;-----------------------------------------------------   
;-----------------------------------------------------   
;-----------------------------------------------------   
;-----------------------------------------------------   
;pipeline swallows first 5 images, which re-numbers everything.
for ca = 1 , n_elements(brightaor) - 1 do begin

   command1 = ' find '+dir_name+'/bcd/'+brightaor(ca)+' -name "*.fits" > ' + dir_name + '/bcd/'+brightaor(ca)+'/files_test.list'
   command2 = 'grep IRAC.1 < '+dir_name+'/bcd/'+brightaor(ca)+'/files_test.list  | grep -v IRAC.2 | grep bcd_fp.fits > '+dir_name+'/bcd/'+brightaor(ca)+'/ch1_bcd_test.list'
   command3 = 'grep IRAC.2 < '+dir_name+'/bcd/'+brightaor(ca)+'/files_test.list  | grep bcd_fp.fits > '+dir_name+'/bcd/'+brightaor(ca)+'/ch2_bcd_test.list'

   a = [command1, command2, command3]
   for i = 0, n_elements(a) -1 do spawn, a(i)

   readcol, dir_name+'/bcd/'+brightaor(ca)+'/ch1_bcd_test.list', bcdname_ch1, format="A",/silent
   readcol, dir_name+'/bcd/'+brightaor(ca)+'/ch2_bcd_test.list', bcdname_ch2, format="A",/silent
   
  
;-----------------------------------------------------   
;forchannel 1
   k = 0
   bcdname = bcdname_ch1

   flux_star1_pos1 = fltarr(20)
   flux_star1_pos2 = fltarr(12)
   flux_star2_pos1 = fltarr(20)
   flux_star2_pos2 = fltarr(12)
   
   
   imagenum = 0
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos1(0) = flux(0)
   flux_star2_pos1(0) = flux(1)
   
   getrot, header, rot, cdelt              ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                     ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0)    ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])           ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)               ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1)    ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])           ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)               ;Within 5 arc seconds
   
   for imagenum = 1, 18 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux
 
;aper is being a shit and won't work for the positions without
;a star on them

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos1(imagenum ) = flux1
      flux_star2_pos1(imagenum ) = flux2
   endfor
   ;------------------

   imagenum = 8
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos2(0) = flux(0)
   flux_star2_pos2(0) = flux(1)
   
   getrot, header, rot, cdelt                  ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                         ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0) ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])             ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)                 ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1) ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])             ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)                 ;Within 5 arc seconds

   for imagenum = 9, 18 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos2(imagenum -8) = flux1
      flux_star2_pos2(imagenum -8) = flux2
   endfor

plot, findgen(20), flux_star1_pos1 / max(flux_star1_pos1), /ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch1 star1 12s'
oplot, findgen(12), flux_star1_pos2 / max(flux_star1_pos2), thick=3, color=bluecolor

plot, findgen(20), flux_star2_pos1/ max(flux_star2_pos1),/ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch1 star2 12s'
oplot, findgen(12), flux_star2_pos2/ max(flux_star2_pos2), thick=3, color=bluecolor

;-----------------------------------------------------   
;forchannel 2
   k = 1
   bcdname = bcdname_ch2

   flux_star1_pos1 = fltarr(24)
   flux_star1_pos2 = fltarr(16)
   flux_star1_pos3 = fltarr(8)
   flux_star2_pos1 = fltarr(24)
   flux_star2_pos2 = fltarr(16)
   flux_star2_pos3 = fltarr(8)
   
   
   imagenum = 4
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos1(0) = flux(0)
   flux_star2_pos1(0) = flux(1)
   
   getrot, header, rot, cdelt              ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                     ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0)    ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])           ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)               ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1)    ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])           ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)               ;Within 5 arc seconds
   
   for imagenum = 5, 18 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
      
;aper is being a shit and won't work for the positions without
;a star on them

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos1(imagenum -4) = flux1
      flux_star2_pos1(imagenum -4) = flux2
   endfor
   ;------------------

   imagenum = 12
   fits_read, bcdname(imagenum), data, header
   adxy, header, ra, dec, xcen, ycen
   
   aper,data, xcen, ycen, flux , fluxerr, sky, $
        skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
   flux_star1_pos2(0) = flux(0)
   flux_star2_pos2(0) = flux(1)
   
   getrot, header, rot, cdelt                  ;CDELT gives plate scale deg/pixel
   cdelt = cdelt*3600.                         ;Convert to arc sec/pixel
   dist_circle, circle0, 256, xcen(0), ycen(0) ;Create a distance circle image
   circle0 = circle0*abs(cdelt[0])             ;Distances now given in arcseconds
   good0 = where(circle0 LT 6)                 ;Within 5 arc secons

   dist_circle, circle1, 256, xcen(1), ycen(1) ;Create a distance circle image
   circle1 = circle1*abs(cdelt[0])             ;Distances now given in arcseconds
   good1 = where(circle1 LT 6)                 ;Within 5 arc seconds

   for imagenum = 13, 18 do begin
      fits_read, bcdname(imagenum), data, header
      aper,data, xcen, ycen, flux , fluxerr, sky, skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

      flux1 =  total( data[good0] )  - sky*n_elements(data[good0])
      flux2 = total( data[good1] )  - sky*n_elements(data[good1])
      for i =0, n_elements(flux1) -1 do if flux1(i) lt 0 then flux1(i) =1E-6
      for i =0, n_elements(flux2) -1 do if flux2(i) lt 0 then flux2(i) =1E-6

      flux_star1_pos2(imagenum -12) = flux1
      flux_star2_pos2(imagenum -12) = flux2
   endfor

  
plot, findgen(24), flux_star1_pos1 / max(flux_star1_pos1), /ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch2 star1 12s'
oplot, findgen(16), flux_star1_pos2 / max(flux_star1_pos2), thick=3, color=bluecolor
;oplot, findgen(8), flux_star1_pos3 / max(flux_star1_pos3), thick=3, color=redcolor

plot, findgen(24), flux_star2_pos1/ max(flux_star2_pos1),/ylog, yrange=[0.01, 1], thick =3, charthick=3, xthick=3,ythick=3, xtitle = 'frame number', ytitle = 'fractional flux', title = 'ch2 star2 12s'
oplot, findgen(16), flux_star2_pos2/ max(flux_star2_pos2), thick=3, color=bluecolor
;oplot, findgen(8), flux_star2_pos3/ max(flux_star2_pos3), thick=3, color=redcolor

;-----------------------------------------------------   

   
endfor


ps_close, /noprint,/noid

end
