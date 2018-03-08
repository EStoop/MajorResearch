pro make_staring_dark
 
 chnum =2

 ;need to prepare flats of the right dimension in both channels [32,32,64]
  caldir_ch1 = '/Users/jkrick/irac_warm/darks/staring/r51660288/ch1/cal/'
  fits_read, strcompress(caldir_ch1+ 'irac_b1_sa_superskyflat_100426.fits',/remove_all), flatdata_ch1, flatheader
  flat64_ch1 = fltarr(32,32,64)
  flatsingle = flatdata_ch1[*,*,0]
  for f = 0, 63 do flat64_ch1[*,*,f] = flatsingle
  caldir_ch2 = '/Users/jkrick/irac_warm/darks/staring/r51677440/ch2/cal/'
  fits_read, strcompress(caldir_ch2+ 'irac_b2_sa_superskyflat_100426.fits',/remove_all), flatdata_ch2, flatheader
  flat64_ch2 = fltarr(32,32,64)
  flatsingle = flatdata_ch2[*,*,0]
  for f = 0, 63 do flat64_ch2[*,*,f] = flatsingle
 
  expname = ['s2s', 's0p4s', 's0p1s', 's0p02s']

  for e = 0,  n_elements(expname) -1 do begin
     case expname(e) OF
        's0p02s': exptime = 0.02
        's0p1s': exptime = 0.1
        's0p4s': exptime = 0.4
        's2s': exptime = 2.0
     endcase
     
     dirname = '/Users/jkrick/irac_warm/darks/staring/'
     cd, dirname
     if chnum eq 1 then command1 =  ' ls r*/ch1*/bcd/*bcd.fits >  darklist.txt'  ; includes all frametimes
     if chnum eq 2 then command1 =  ' ls r*/ch2*/bcd/*bcd.fits >  darklist.txt'  ; includes all frametimes
     spawn, command1
     readcol,'darklist.txt',fitsname, format = 'A', /silent
     print, 'n darks', n_elements(fitsname)
     
      
;ugh 4 dimensions
     bigim = fltarr(32, 32, 64, n_elements(fitsname))
     count = 0
     
     for i =0, n_elements(fitsname) - 1 do begin ;
        header = headfits(fitsname(i))           ;
        ch = sxpar(header, 'CHNLNUM')
        framtime = sxpar(header, 'FRAMTIME')
        naxis = sxpar(header, 'NAXIS')
        fluxconv = sxpar(header, 'FLUXCONV')

                                ;make sure I really pulled the correct darks
        if ch eq chnum and naxis eq 3 and framtime eq exptime then begin
           fits_read, fitsname(i), data, header
           data = data / fluxconv ; now in DN/s
           data = data* exptime   ; now in DN
                                ; flip the image
           data = reverse(data, 2)
           
                                ;undo the flat removal
           if chnum eq 1 then data = data * flat64_ch1
           if chnum eq 2 then data = data * flat64_ch2

                                ;remove the dark that was already used in the image
           darkname = sxpar(header, 'SKDKRKEY')
           darkepid = sxpar(header, 'SDRKEPID')
           framedelay = sxpar(header, 'FRAMEDLY')
           aorkey = sxpar(header, 'AORKEY')
           darkname = strcompress('r'+string(aorkey) + '/ch' + string(chnum)+'/cal/SPITZER_I'+string(chnum)+'_'+string(darkname)+ '_00*_C'+string(darkepid)+'_sdark.fits',/remove_all)
           fits_read,  darkname, darkdata, darkheader
           data = data + darkdata
           
           ;now fill the bigim array
           bigim(0, 0, 0, count) = data
           count = count + 1
           if i eq n_elements(fitsname) - 1 then savedata = data
        endif
     endfor
     print, 'count', count
     if count ne n_elements(fitsname) then bigim= bigim[*,*,*,0:count-1]
     
;now make a median
     staringdark = median(bigim, dimension = 4)
     fits_write, strcompress(dirname+ 'staringdark_'+expname(e)+'_ch'+ string(chnum) +'.fits',/remove_all), staringdark, header ; just use the last header
     
  endfor   ; for each expname


end
