pro mean_staring_dark
  
;what is the mean level as a function of frame number and amplifier and frame delay

;  aorname = 'r40151040'        ; this one has a dark change in the middle of it.
  ;40151552 has a star in it, 
  aorname = ['/Users/jkrick/irac_warm/hd189733/r40152064/' ,'/Users/jkrick/irac_warm/hd189733/r40151040/'  , '/Users/jkrick/irac_warm/hd189733/r40150016/', '/Users/jkrick/irac_warm/hd189733/r40150528/', '/Users/jkrick/irac_warm/hd189733/r40151552/', '/Users/jkrick/irac_warm/HD209458/r44201728/', '/Users/jkrick/irac_warm/HD209458/r44201216/', '/Users/jkrick/irac_warm/HD209458/r44201472/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42051584/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42506496/']
  
  chname =[  '2',  '2', '2', '2', '2', '1',  '1', '1', '1', '1']
  
  a = findgen(32)
  b = a mod 4
  amp0 = where(b eq 0)
  amp1 = where(b eq 1)
  amp2 = where(b eq 2)
  amp3 = where(b eq 3)
  
                                ;for each aor I need to carry around
 ; fitsname, framenumber(0: n_elements(fitsname)), imagenum[0:63,0:63,...)], mean & stddev in each amplifier[0:63], framedelay[0:n_elements(fitsname)]
  biashash = hash()
 
  for a = 0, 0 do begin; n_elements(aorname) - 1 do begin
     print, 'working on aor', aorname(a)
     skyflatname = strcompress( aorname(a) + 'ch' + chname(a) + '/cal/*superskyflat*.fits',/remove_all)
;     print, 'skyflatname', skyflatname
     fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
     flat64 = fltarr(32,32,64)
     flatsingle = flatdata[*,*,0]
     for f = 0, 63 do flat64[*,*,f] = flatsingle
     
    
     CD, aorname(a)             ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname(a)+"/bcd -name 'SPITZER*_bcd.fits' > "+aorname(a)+'bcdlist.txt')
     spawn, command
     readcol,strcompress(aorname(a) +'bcdlist.txt'),fitsname, format = 'A', /silent
     
     ;arrays over all subframes in the AOR
     meanarr_0 = fltarr(n_elements(fitsname)) ; total number of pixels in a single amplifier over whole AOR
     meanarr_1 = meanarr_0
     meanarr_2 = meanarr_0
     meanarr_3 = meanarr_0
     sigmaarr_0 = meanarr_0
     
     for i = 0,  n_elements(fitsname) - 1,5 do begin
        fits_read, fitsname(i), data, header
        
                                ;back out the flux conversion
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        data = data / fluxconv  ; now in DN/s
        data = data* exptime  ; now in DN
                                ; flip the image
        data = reverse(data, 2)
        
                                ;remove the flat
        data = data * flat64
        
        darkname = sxpar(header, 'SKDKRKEY')
        darkepid = sxpar(header, 'SDRKEPID')
        framedelay = sxpar(header, 'FRAMEDLY')
        aorkey = sxpar(header, 'AORKEY')

                                ;remove the dark that was already used in the image
        fits_read, strcompress(aorname(a) + 'ch' + chname(a)+'/cal/SPITZER_I'+chname(a)+'_'+string(darkname)+$
                               '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), darkdata, darkheader
        data = data + darkdata
        
        ;just one amplifier for now
        meanclip, data[amp0,*,*], m0, s0
        meanarr_0(i) =  m0
        sigmaarr_0(i) = s0
     endfor ; for each fitsname

    
     ;plotting
     x = indgen(n_elements(meanarr_0))
     p = errorplot(x, meanarr_0, sigmaarr_0,  '1s', xtitle = 'Frame number', ytitle = 'mean level', sym_size = 0.5, yrange = [0,20], xrange = [0,650])

    endfor                        ; for each aorname
  
end
