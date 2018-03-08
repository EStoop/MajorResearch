pro pmap_darkcorr, frametime, chname
;this code backs the dark correction out of the snapshot observations,
;and then puts back the same dark for all snaps.

;assumes subarray - aka naxis = 3
  

  basedir = '/Users/jkrick/irac_warm/calstars/pmap_star_ch2/'
  aorname = 'r46037248'
  
 
; the superdark
  if frametime eq 2 then frametimename = 's2s'
  if frametime eq 0.1 then frametimename = 's0p1s'
  fits_read,strcompress( '/Users/jkrick/external/irac_warm/darks/superdarks/superdark_'+frametimename +'.fits',/remove_all), superdark, superheader

                      
 
  prefluxarr = fltarr(n_elements(aorname) * 500* 64)  ;just guesses at size for now
  superfluxarr = prefluxarr
  c = 0L

  for a = 0,   n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = basedir + string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dir+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name 'SPITZER*_bunc.fits' > "+dir+'bunclist.txt')
     spawn, command2
   
     readcol,strcompress(dir +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dir +'bunclist.txt'),buncname, format = 'A', /silent
     print,'n_elements(fitsname)', n_elements(fitsname)
     
       ;need to read in the flat and make this [32,32,64]
     skyflatname = strcompress( 'ch' + chname + '/cal/*superskyflat*.fits',/remove_all)
     fits_read, skyflatname, flatdata, flatheader
     flat64 = fltarr(32,32,64)
     flatsingle = flatdata[*,*,0]
     for f = 0, 63 do flat64[*,*,f] = flatsingle
 

     for i = 0,   n_elements(fitsname) - 1 do begin ; for each image

                                ;read in the data
        fits_read, fitsname(i), data, header
        fits_read, buncname(i), unc, uncheader

        ra_ref = sxpar(header, 'RA_RQST')
        dec_ref = sxpar(header, 'DEC_RQST')

        ; do some baseline aperture photometry
        get_centroids_for_calstar_jk,data, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        prefluxarr(c) = f[*,1]      
        
                                ;back out the flux conversion
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        data = data / fluxconv  ; now in DN/s
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
        
                                ;write out the new fits file
        outfilename =strmid(fitsname(i),0, 8) + 'sdcorr' + strmid(fitsname(i), 8)
        fits_write, outfilename, data, header
        

        ;do some aperture photometry after changing the dark
        get_centroids_for_calstar_jk,data, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        superfluxarr(c) = f[*,1]      
        c = c + 63

     endfor                     ; for each image
  endfor                        ; for each AOR

  prefluxarr = prefluxarr[0:c-63]
  superfluxarr = superfluxarr[0:c-63]
  deltaflux = prefluxarr - superfluxarr
  binsize = 0.00005
  plothist, deltaflux, prexhist, preyhist, /noplot, bin = binsize
  h = plot(prexhist, preyhist, thick = 2, xtitle = 'delta flux', ytitle = 'Number', title = aorname(0), $
           color = 'Navy', name = 'campaign dark', xrange = [median(deltaflux) - 20.*binsize, median(deltaflux) + 20.*binsize])
  
;  plothist, superfluxarr, superxhist, superyhist, /noplot, bin = binsize
;  h2 =  plot(superxhist, superyhist, thick = 2, /overplot, color = 'red', name = 'superdark')

;  l = legend(target = [h, h2], position = [median(prefluxarr), 1000], /data, /auto_text_color)

end
