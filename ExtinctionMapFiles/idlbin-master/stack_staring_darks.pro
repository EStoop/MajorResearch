pro stack_staring_darks, ch


;  aorname = 'r40151040'        ; this one has a dark change in the middle of it.
  restore,  '/Users/jkrick/irac_warm/darks/staring/staring_darks.sav'
  aorname = ['/Users/jkrick/irac_warm/hd189733/r40152064/' ,'/Users/jkrick/irac_warm/hd189733/r40151040/'  , '/Users/jkrick/irac_warm/hd189733/r40150016/', '/Users/jkrick/irac_warm/hd189733/r40150528/', '/Users/jkrick/irac_warm/hd189733/r40151552/', '/Users/jkrick/irac_warm/HD209458/r44201728/', '/Users/jkrick/irac_warm/HD209458/r44201216/', '/Users/jkrick/irac_warm/HD209458/r44201472/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42051584/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42506496/']
  
  chname =[  '2',  '2', '2', '2', '2', '1',  '1', '1', '1', '1']
    
  if ch eq 2 then begin
     startaor = 0
     stopaor = 4      
     meanyrange = [9.5,11.5]
  endif
  if ch eq 1 then begin
     startaor = 5; n_elements(aorname) - 1
     stopaor = n_elements(aorname) - 1
     meanyrange = [116, 130]
  endif
  
  for a = startaor,  stopaor do begin
     print, 'working on aor', aorname(a)
     skyflatname = strcompress( aorname(a) + 'ch' + chname(a) + '/cal/*superskyflat*.fits',/remove_all)
     fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
     flat64 = fltarr(32,32,64)
     flatsingle = flatdata[*,*,0]
     for f = 0, 63 do flat64[*,*,f] = flatsingle
    
     CD, aorname(a)             ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname(a)+"/bcd -name 'SPITZER*_bcd.fits' > "+aorname(a)+'bcdlist.txt')
     spawn, command
     readcol,strcompress(aorname(a) +'bcdlist.txt'),fitsname, format = 'A', /silent
     
     ;want this to run from 0 to 63 for n_elements(fitsname)
;     for si = 0, n_elements(fitsname)*63, 64 do subimagenum(si) = findgen(64)
;     c = 0L
     
     nstack = 50 ; n_elements(fitsname) - 1

     ; make a stack with images 0 - nstack, and 550 + nstack
     medimearly = make_stack(ch, fitsname, flat64, 0, nstack)
     medimlate = make_stack(ch, fitsname, flat64, 550, nstack)
     print, 'medimearly', medimearly[30, 14, 14]
     print, 'medimlate', medimlate[30, 14, 14]

     ;subtract them, see what comes out
     diffim = medimearly - medimlate
     fits_write,  strcompress( '/Users/jkrick/irac_warm/darks/staring/diffim_'+string(biashash[aorname(a),'aorkey'])+'.fits',/remove_all), diffim, header ; just use the last header
     
     print, 'stddev diffim', stddev(diffim)
     mmm, diffim, skymod, sigma, skew,/silent   ; should do a better clipping job of outliers
     print, 'mmm sigma', sigma

     ;try aperture photometry of the center pixel
;     for j = 0, 63 do begin
        
;        aper, diffim[*,*,j], 15.5, 15.5, flux, fluxerr, back, backerr, 1.0, 3, [7,10], $
;              9E8, /FLUX, /EXACT, /SILENT, /NAN, /MEANBACK
;        print, 'aper', flux, fluxerr, back, backerr
 ;    endfor

  endfor   ; for each aorname


  
end


function make_stack, ch, fitsname, flat64,  nstart, nstack
    ;set up the big array for stacking
  bigim = fltarr(32, 32, 64, nstack )
  count = 0
  
  for i = nstart,  nstart + nstack - 1 do begin
     fits_read, fitsname(i), data, header
     
                                ;back out the flux conversion
     fluxconv = sxpar(header, 'FLUXCONV')
     exptime = sxpar(header, 'EXPTIME')
     data = data / fluxconv     ; now in DN/s
     data = data* exptime       ; now in DN
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
     
                                ;stack-em
     bigim(0, 0, 0, count) = data
     
;        help, bigim
;        help, data
                                ;print, 'bigim as its built', bigim[14,14,14, *]
     count = count + 1
  endfor                        ; for each fits image in desired stack
  
                                ;now median the whole stack together
  medim = median(bigim, dimension = 4)
  
  return, medim
end
