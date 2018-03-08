pro singlepix_staring_darks, ch
;  aorname = 'r40151040'        ; this one has a dark change in the middle of it.
  restore,  '/Users/jkrick/irac_warm/darks/staring/staring_darks.sav'
  aorname = ['/Users/jkrick/irac_warm/hd189733/r40152064/' ,'/Users/jkrick/irac_warm/hd189733/r40151040/'  , '/Users/jkrick/irac_warm/hd189733/r40150016/', '/Users/jkrick/irac_warm/hd189733/r40150528/', '/Users/jkrick/irac_warm/hd189733/r40151552/', '/Users/jkrick/irac_warm/HD209458/r44201728/', '/Users/jkrick/irac_warm/HD209458/r44201216/', '/Users/jkrick/irac_warm/HD209458/r44201472/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42051584/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42506496/']
  colorarr = ['black', 'red', 'blue', 'cyan', 'purple', 'black', 'red', 'blue', 'cyan', 'purple']
  chname =[  '2',  '2', '2', '2', '2', '1',  '1', '1', '1', '1']
    
  if ch eq 2 then begin
     startaor = 0
     stopaor = 4      
     meanyrange = [9.5,11.5]
     xr = [0,650]
     yr = [0.6, 1.4]
  endif
  if ch eq 1 then begin
     startaor = 5               ; n_elements(aorname) - 1
     stopaor = n_elements(aorname) - 1
     meanyrange = [116, 130]
     xr = [0,1300]
     yr = [0.96, 1.03]
  endif
  chname = strcompress(string(ch),/remove_all)

  for a = startaor,  stopaor do begin
     print, 'working on aor', aorname(a)
     skyflatname = strcompress( aorname(a) + 'ch' + chname + '/cal/*superskyflat*.fits',/remove_all)
     fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
     flat64 = fltarr(32,32,64)
     flatsingle = flatdata[*,*,0]
     for f = 0, 63 do flat64[*,*,f] = flatsingle
    
     CD, aorname(a)             ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+aorname(a)+'bcdlist.txt')
     spawn, command
     readcol,strcompress(aorname(a) +'bcdlist.txt'),fitsname, format = 'A', /silent
     
     pix1arr = fltarr(n_elements(fitsname))
     for i = 0,  n_elements(fitsname)- 1 do begin
        fits_read, fitsname(i), data, header
        
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
        darkname = sxpar(header, 'SKDKRKEY')
        darkepid = sxpar(header, 'SDRKEPID')
        framedelay = sxpar(header, 'FRAMEDLY')
        aorkey = sxpar(header, 'AORKEY')
        chname = strcompress(string(ch),/remove_all)
        fits_read, strcompress( 'ch' + chname+'/cal/SPITZER_I'+chname+'_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), $
                   darkdata, darkheader
        data = data + darkdata
        
        
                                ;follow one pixel on one subframe
        pix1arr(i) = data[9,12,10]
     endfor  ; each fits image
     
     ;normalize by the mean 
     normfactor = mean(pix1arr,/nan)
     pix1arr = pix1arr / normfactor

     ;boxcar average over N pixels
     nsmooth = 10
     smootharr = smooth(pix1arr, nsmooth,/NAN,/EDGE_TRUNCATE)
     p = plot(indgen(n_elements(smootharr)), smootharr, xtitle = 'image number', ytitle = 'single pixel value', $
              title ='Smoothing by '+string(nsmooth) + '  ch'+ string(ch), color = colorarr(a), xrange = xr, yrange = yr, /overplot)

  endfor ; each aor
end
