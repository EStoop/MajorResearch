pro run_latent, dir_name, aor_bright

;ps_open, filename='/Users/jkrick/iwic/iwic_recovery12/IRAC018400/latent.ps',/portrait,/square,/color

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

for ca = 0 , n_elements(brightaor) - 1 do begin

   command1 = ' find '+dir_name+'/bcd/'+brightaor(ca)+' -name "*.fits" > ' + dir_name + '/bcd/'+brightaor(ca)+'/files_test.list'
   command2 = 'grep IRAC.1 < '+dir_name+'/bcd/'+brightaor(ca)+'/files_test.list  | grep bcd_fp.fits > '+dir_name+'/bcd/'+brightaor(ca)+'/ch1_bcd_test.list' ;grep -v IRAC.2 | 
   command3 = 'grep IRAC.2 < '+dir_name+'/bcd/'+brightaor(ca)+'/files_test.list  | grep bcd_fp.fits > '+dir_name+'/bcd/'+brightaor(ca)+'/ch2_bcd_test.list'

   a = [command1, command2, command3]
   for i = 0, n_elements(a) -1 do spawn, a(i)

   readcol, dir_name+'/bcd/'+brightaor(ca)+'/ch1_bcd_test.list', bcdname_ch1, format="A",/silent
   readcol, dir_name+'/bcd/'+brightaor(ca)+'/ch2_bcd_test.list', bcdname_ch2, format="A",/silent

  
   nimages = n_elements(bcdname_ch1)
   
;two stars to keep track of
   xcenarr_1 = (fltarr(nimages))
   ycenarr_1 = fltarr(nimages)
   xcenarr_2 = (fltarr(nimages))
   ycenarr_2 = fltarr(nimages)
   
;for each channel
   for k = 0, 1 do begin
      if k eq 0 then bcdname = bcdname_ch1
      if k eq 1 then bcdname = bcdname_ch2
;find the locations of the stars on all frames
      for i = 0, nimages - 1 do begin

         header = headfits(bcdname(i)) ;read in the header of the image on which we are working
;         print, 'frametime',  fxpar(header, 'framtime')

         adxy, header, ra, dec, xcen, ycen
;         print, 'xcen, ycen', xcen, ycen
                                ;keep track of where the two stars
                                ;are, only keep if they are on the image
         xcenarr_1(i) = xcen(0)
         ycenarr_1(i) = ycen(0)
         xcenarr_2(i) = xcen(1)
         ycenarr_2(i) = ycen(1)
         
      endfor

 
      fluxarr_1 = fltarr(nimages,nimages)
      fluxarr_2 = fltarr(nimages,nimages)
      
;measure the flux at each of those positions
      for i = 0, nimages - 1 do begin
;      print, 'working on ', bcdname(i)
         
         fits_read, bcdname(i), data, header
         
         ;won't work for subarray, ignore those
         naxis1 = fxpar(header, 'NAXIS1')
         if naxis1 gt 50 then begin
         ;star 1. on one image measure the flux at all possible star one positions
            aper,data, xcenarr_1, ycenarr_1, flux , fluxerr, sky, $
                 skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
            flux(where(finite(flux) lt 1)) = 0
            fluxarr_1(i,*) = flux

            ;star 2
            aper,data, xcenarr_2, ycenarr_2, flux2 , fluxerr2, sky2, $
                 skyerr2, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
            nanflux = where(finite(flux2) lt 1)
            flux2(where(finite(flux2) lt 1)) = 0

            fluxarr_2(i,*) = flux2
         endif
         
      endfor
      
      
   ;plot flux as a function of frame number (time)
   ;star 1
      print, 'test', fluxarr_1(0,*)
      plot, findgen(nimages), fluxarr_1(0,*) / max(fluxarr_1(0,*)), yrange=[1E-4,1], thick = 3, charthick = 3, /ylog,$
            xthick = 3, ythick = 3, title =  strcompress('Ch' + string( k+1) + ' STAR 1'),$
            xtitle = 'Frame Number', ytitle = 'Fractional Flux '
      oplot, findgen(nimages), fluxarr_1(0,*)/ max(fluxarr_1(0,*)), psym = 2, thick = 3
      
      oplot, findgen(nimages), fluxarr_1(1,*) / max(fluxarr_1(1,*)), thick = 3, color = colorname[1]
      oplot, findgen(nimages), fluxarr_1(1,*) / max(fluxarr_1(1,*)), thick = 3, color = colorname[1], psym = 2
      
      for i =2, nimages - 2 do begin
         if i mod 2 eq 0 then newi = i
         if i  mod 2 eq 1 then newi = i -1
         oplot, findgen(nimages-2), fluxarr_1(i,newi:*) / max(fluxarr_1(i,*)), thick = 3, color = colorname[i]
         oplot, findgen(nimages-2), fluxarr_1(i,newi:*) / max(fluxarr_1(i,*)), thick = 3, color = colorname[i], psym = 2
      endfor
      
  ;print out the fluxes at the star to see if photometry is stable.
      y = [.013,.02, .03,.05,.08, .13, .2,.3,.5,.8]
;      for i = 0, nimages - 2 do xyouts, 6, y(i), max(fluxarr_1(i,*))  

   ;star 2
      print, 'test', fluxarr_2(0,*)
      
      plot, findgen(nimages), fluxarr_2(0,*) / max(fluxarr_2(0,*)), yrange= [1E-4,1], thick = 3, charthick = 3, $
            xthick = 3, ythick = 3, title =  strcompress('Ch' + string( k+1) + ' STAR 2'), xtitle = 'frame number', ytitle = 'Fractional Flux ', /ylog
      oplot, findgen(nimages), fluxarr_2(0,*) / max(fluxarr_2(0,*)), psym = 2, thick = 3
      
      oplot, findgen(nimages), fluxarr_2(1,*) / max(fluxarr_2(1,*)), thick = 3, color = colorname[1]
      oplot, findgen(nimages), fluxarr_2(1,*) / max(fluxarr_2(1,*)), thick = 3, color = colorname[1], psym = 2
      

      for i =2, nimages - 2 do begin ;2
         if i mod 2 eq 0 then newi = i
         if i  mod 2 eq 1 then newi = i -1
         oplot, findgen(nimages-2), fluxarr_2(i,newi:*) / max(fluxarr_2(i,*)), thick = 3, color = colorname[i]
         oplot, findgen(nimages-2), fluxarr_2(i,newi:*) / max(fluxarr_2(i,*)), thick = 3, color = colorname[i], psym = 2
      endfor
      
;;      for i = 0, nimages - 2 do xyouts, 6, y(i), max(fluxarr_2(i,*))

   endfor


;----------------

;display them in ds9
   
   for j = 0,  1 do begin       ;for each channel
      if j eq 0 then bcdname = bcdname_ch1
      if j eq 1 then bcdname = bcdname_ch2
      command = 'ds9 '
      for i = 0, n_elements(bcdname) -1 do  begin
         command = command + ' ' + bcdname(i) + ' -zscale '
                                ;command = command + ' -single'
      endfor

      print, 'A ds9 image is opening for channel ' + string(j+ 1)
      print, 'Look at the bright sources images for artifacts or latents.'
      print, 'Close the ds9 window to continue.'
      latents = 'no'
   while (latents eq 'no') do begin
      spawn, command
      read, latents, prompt = 'Are you done looking for artifacts in the ds9 window? (yes or no) '
   endwhile


   endfor
endfor

;ps_close, /noprint,/noid

end
