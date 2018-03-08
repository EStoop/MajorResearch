pro phot_pmap, planetname, apradius,chname, breatheap = breatheap
;do photometry on  IRAC pmap dataset
 t1 = systime(1)

;convert aperture radius in pixels into what get_centroids_for_calstar_jk uses 
case apradius of
   1.5: apval = 0
   1.75: apval = 1
   2.0: apval = 2
   2.25:  apval = 1
   2.5: apval = 2
   2.75: apval = 5
   3.25: apval = 7
   3: apval = 6
   Else: apval = 2              ; if no decision, then choose an apradius = 2.5 pixels
endcase


;run code to read in all the planet parameters
;if chname eq '2'  then aorname= aorname_ch2 else aorname = aorname_ch1
aorname = ['r38933504', 'r38933760', 'r38934016', 'r38934272', 'r38934528', 'r38934784', 'r38935040', 'r38935296', 'r38936064', 'r38936320', 'r38935552', 'r38936576', 'r38935808', 'r38936832', 'r39112960', 'r39113472', 'r39113216', 'r39113728', 'r39113984', 'r39115264', 'r39114240', 'r39114496', 'r39114752', 'r39115008', 'r39115520', 'r39115776', 'r39116032', 'r39116288']
ra_ref = 13  ;XXX
dec_ref = 13

dirname = '/Users/jkrick/irac_warm/pmap/'
planethash = hash()


startaor =0
stopaor = n_elements(aorname) - 1
for a =startaor, stopaor do begin
   print, 'working on ',aorname(a)
   dir = dirname+ string(aorname(a) ) 
   CD, dir                      ; change directories to the correct AOR directory
   command  = strcompress( 'find ch'+chname+"/bcd -name '*bcd.fits' > "+dirname+'bcdlist.txt')
   print, 'command', command
   spawn, command
   command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
   spawn, command2
   
   readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
   readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
   
   print,'n_elements(fitsname)', n_elements(fitsname)
;     aparr = dblarr(n_elements(fitsname))  ;keep the aperture sizes used
  

       ;need to read in the flat and make this [32,32,64]
   skyflatname = strcompress( 'ch' + chname + '/cal/*superskyflat*.fits',/remove_all)
   fits_read, skyflatname, flatdata, flatheader
   flat64 = fltarr(32,32,64)
   flatsingle = flatdata[*,*,0]
   for f = 0, 63 do flat64[*,*,f] = flatsingle


   for i =0.D,  n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
 ;      print, 'working on ', fitsname(i)         
      header = headfits(fitsname(i)) ;
      sclk_obs= sxpar(header, 'SCLK_OBS')
      frametime = sxpar(header, 'FRAMTIME')
      bmjd_obs = sxpar(header, 'BMJD_OBS')
       ;utcs_obs = sxpar(header, 'UTCS_OBS')
      ch = sxpar(header, 'CHNLNUM')
      ronoise = sxpar(header, 'RONOISE') ; these are zeros
      gain = sxpar(header, 'GAIN')
      fluxconv = sxpar(header, 'FLUXCONV')
      exptime = sxpar(header, 'EXPTIME')
      aintbeg = sxpar(header, 'AINTBEG')
      atimeend = sxpar(header, 'ATIMEEND')
      naxis = sxpar(header, 'NAXIS')
      if ch eq '2' and frametime eq 2 then ronoise = 12.1
      if i eq 0 then print, 'ronoise', ronoise, gain, fluxconv, exptime
      if i eq 0 then sclk_0 = sclk_obs

      ;setup arrays to hold the photometry
      if i eq 0 and naxis eq 3 then begin
         xarr = fltarr(63*(n_elements(fitsname)))
         bmjd = dblarr(63*(n_elements(fitsname)))
      endif
      if i eq 0 and naxis ne 3 then begin
         xarr = fltarr(n_elements(fitsname))
         bmjd = dblarr(n_elements(fitsname))
       endif
      yarr = xarr
      fluxarr = xarr
      fluxerrarr = xarr
      timearr = xarr
      backarr = xarr
      backerrarr = xarr
      nparr = xarr
      npcentroidsarr = xarr
      
      
      if naxis eq 3 then begin
         deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
         nt = dindgen(64)
         sclkarr = sclk_obs  + (deltatime*nt)/60./60./24.D ; 0.5*frametime + frametime*nt
         bmjdarr= bmjd_obs + (deltatime*nt)/60./60./24.D   ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
      endif else begin          ; full array, so don't need to expand out the times
         sclkarr = sclk_obs
         bmjdarr = bmjd_obs
      endelse

      ;read in the files
      fits_read, fitsname(i), im, h
      fits_read, buncname(i), unc, hunc
      test = sxpar(h, 'RAWFILE')
      
      ;XXXXXX
                                ;need to make a superdark
                                ;need to figure out which AORs have which channel, and the ra and dec
      ; wait these are full array?

                                ;back out the flux conversion
      data = im / fluxconv    ; now in DN/s
        data = data* exptime    ; now in DN
        
                                ; flip the image
        data = reverse(data, 2)
        
                                ;remove the flat
        data = data * flat64
        
                                ;remove the dark that was already used in the image
        command  = strcompress( 'find ch'+chname+"/cal -name '*_sdark.fits' > "+dir+'dark.txt')
        spawn, command
        readcol,strcompress(dir +'dark.txt'),darkname, format = 'A', /silent
        fits_read, darkname, dark, darkheader
        
        data = data + dark
       
                                ;--------------------------
        
                                ;put the same dark back into every frame
        data = data - superdark
        
                                ;divde back out the flat
        data = data / flat64
        
                                ;flip the image
        data = reverse(data, 2)
        
                                ;put the fluxconv back in
        data = data / exptime   ; DN/s
        data = data * fluxconv  ;MJy/sr


        ;run the centroiding and photometry
      get_centroids_for_calstar_jk,data, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm, /WARM
      x_center = temporary(x3)
      y_center = temporary(y3)
     ;choose the requested pixel aperture
      abcdflux = f[*,apval]      
      fs = fs[*,apval]
     ; 3-7 pixel background
      back = b[*,0]
      backerr = bs[*,0]
      npcentroids = np   ;keep track of np from get_centroids

      ; if changing the apertures then use this to calculate photometry
      if keyword_set(breatheap) then begin
         abcdflux = betap(im, x_center, y_center, ronoise, gain, exptime, fluxconv,np, chname)
         ;XXXfake these for now
         fs = abcdflux
         back = abcdflux
         backerr = abcdflux
      endif 
      
      


;---------------------------------

      if naxis eq 3 then begin  ; and i eq 0 then begin
         xarr[i*63] = x_center[1:*]
         yarr[i*63] = y_center[1:*]
         fluxarr[i*63] = abcdflux[1:*]
         fluxerrarr[i*63] = fs[1:*]
         corrfluxarr[i*63] = corrflux[1:*]
         corrfluxerrarr[i*63] = corrfluxerr[1:*]
         timearr[i*63] = sclkarr[1:*]        
         bmjd[i*63] = bmjdarr[1:*]
         backarr[i*63] = back[1:*]
         backerrarr[i*63] = backerr[1:*]
         nparr[i*63] = np[1:*]
         npcentroidsarr[i*63] = npcentroids[1:*]
 ;        help, bmjd
      endif 
      if naxis eq 2 then begin; and i eq 0 then begin
         xarr[i] = x_center
         yarr[i]  =  y_center
         fluxarr[i]  =  abcdflux
         fluxerrarr[i]  =  fs
         corrfluxarr[i]  = corrflux
         corrfluxerrarr[i]  =  corrfluxerr
         timearr[i]  = sclkarr
         bmjd[i]  = bmjdarr
         backarr[i]  =  back
         backerrarr[i]  = backerr
         nparr[i]  = np
         npcentroidsarr[i] = npcentroids
      endif

   endfor                       ; for each fits file in the AOR

;--------------------------------
;fill in that hash of hases
;--------------------------------


     keys =['ra', 'dec', 'xcen', 'ycen', 'flux','fluxerr',  'sclktime_0', 'timearr', 'aor', 'bmjdarr', 'bkgd', 'bkgderr', 'npcentroids','exptime']
      values=list(ra_ref,  dec_ref, xarr, yarr, fluxarr, fluxerrarr, sclk_0, timearr, aorname(a), bmjd,  backarr, backerrarr, npcentroidsarr, exptime)
      planethash[aorname(a)] = HASH(keys, values)
 

endfor                          ;for each AOR


if keyword_set(breatheap) then begin
   savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_varap.sav')
endif else begin
   savename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
endelse

save, planethash, filename=savename
print, 'saving planethash', savename
print, 'time check', systime(1) - t1




end





function betap, im, xcen, ycen, ronoise, gain, exptime, fluxconv, np, chname
 ;this function does aperture photometry based on an aperture size that is allowed to vary
  ; as a function of noise pixel
  ; ref: Knutson et al. 2012

  backap = [3.,7.] 
  convfac = gain*exptime/fluxconv
  eim = im * convfac

  varap = sqrt(np)  + 0.2
  ;print, 'testing vartap', varap

  ;XXX add some way of keeping track of varap
  ;don't know how to return that

  abcdflux = fltarr(64,/NOZERO)
  badpix = [-9., 9.] * 1.D8
  pxscal1 = [-1.22334117768332D, -1.21641835430637D, -1.22673962032422D, -1.2244968325831D]
  pxscal1 = abs(pxscal1)
  pxscal2 = [1.22328355209902D, 1.21585676679388D, 1.22298117494211D, 1.21904995758086D]
  pscale2 = pxscal1[long(chname) - 1] * pxscal2[long(chname) - 1]
  scale = pscale2 * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

  for s= 0, 63 do begin
     eslice = eim[*,*,s]
     aper, eslice, xcen[s], ycen[s], xf, xfs, xb, xbs, 1.0, varap[s], backap, $
			      badpix, /FLUX, /EXACT, /SILENT, /NAN, $;
			      READNOISE=ronoise
     f = xf/ convfac
     f = f * scale
     abcdflux[s] =f 
;     print, 'varap, abcdflux', varap[s], abcdflux[s]
  endfor


  return, abcdflux
end

         ;slice up image
;           	for s = 0, 63 do begin
;                   eslice = eim[*, *, s];;;

