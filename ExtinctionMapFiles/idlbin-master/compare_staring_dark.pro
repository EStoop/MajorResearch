pro compare_staring_dark
  
;create median images at each position without the first fits file and without the first 10 subframes frames, and see if there is a pattern change between different frame times.
  dirname = '/Users/jkrick/irac_warm/darks/staring/'
;start w ch1
;0.02s images ; positions 1 - 16
  aor0p02s = ['r51659264',	'r51659520',	'r51659776',	'r51660032',	'r51660288',	'r51660544',	'r51660800',	'r51661056',	'r51661312',	'r51661568',	'r51661824',	'r51662080',	'r51662336',	'r51662592',	'r51662848',	'r51663104']
;0.1s images, positions 1 = 16
  aor0p1s =['r51659008',	'r51663360',	'r51663616',	'r51663872',	'r51664128',	'r51664384',	'r51664640',	'r51664896',	'r51665152',	'r51665408',	'r51665664',	'r51665920',	'r51666176',	'r51666432',	'r51666688',	'r51666944']
;0.4s, positions 1 = 16
  aor0p4s =[	'r51667200',	'r51667456',	'r51667712',	'r51667968',	'r51668224',	'r51668480',	'r51668736',	'r51668992',	'r51669248',	'r51669504',	'r51669760',	'r51670016',	'r51670272',	'r51670528',	'r51670784',	'r51671040']
;2s, positions 1 -2
  aor2s = [	'r51671296',	'r51671552']
  
  allaors = [aor0p02s, aor0p1s, aor0p4s] ;, aor2s]
  chname = '1'
  npositions = 16; which is  = n_elements(aorname)
  expname = ['0p02s','0p1s', '0p4s'];, '2s']
  saveexp = fltarr(32,32,54, npositions, n_elements(expname))
  for e = 0,  n_elements(expname) -1 do begin

     if e eq 0 then aorname = aor0p02s
     if e eq 1 then aorname = aor0p1s
     if e eq 2 then aorname = aor0p4s

     for a = 0, n_elements(aorname) - 1 do begin
        print, 'working on aor: ', aorname(a)
        skyflatname = strcompress( dirname + aorname(a) + '/ch' + chname + '/cal/*superskyflat*.fits',/remove_all)
        fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
        flat64 = fltarr(32,32,64)
        flatsingle = flatdata[*,*,0]
        for f = 0, 63 do flat64[*,*,f] = flatsingle

        CD, dirname + aorname(a) ; change directories to the correct AOR directory
        command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+aorname(a)+'bcdlist.txt')
        spawn, command
        readcol,strcompress(aorname(a) +'bcdlist.txt'),fitsname, format = 'A', /silent
                
                                ;ugh 4 dimensions
        bigim = fltarr(32, 32, 54, n_elements(fitsname))
        
        for i = 1,  n_elements(fitsname) - 1 do begin
           fits_read, fitsname(i), data, header
                                ;back out the flux conversion
           fluxconv = sxpar(header, 'FLUXCONV')
           exptime = sxpar(header, 'EXPTIME')
          data = data / fluxconv ; now in DN/s
           data = data* exptime   ; now in DN
                                ; flip the image
           data = reverse(data, 2)
           
                                ;undo the flat removal
           data = data * flat64
           
                                ;remove the dark that was already used in the image
           darkname = sxpar(header, 'SKDKRKEY')
           darkepid = sxpar(header, 'SDRKEPID')
           framedelay = sxpar(header, 'FRAMEDLY')
           aorkey = sxpar(header, 'AORKEY')
           fits_read, strcompress('ch' + chname+'/cal/SPITZER_I'+chname+'_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), darkdata, darkheader
           data = data + darkdata
           
           ;remove the flat back in !!!???
           data = data / flat64
           
                                ;fill inthe array to make a big median image
           data = data [*,*,10:63]
           bigim(0, 0, 0, i) = data
           
        endfor                  ; for each fitsname
;now make a median
        aordark = median(bigim, dimension = 4)
        saveexp(0,0,0,a, e) = aordark
;        filename = strcompress(dirname+ 'superdark_'+string(a) +'_'+ expname(e)+'.fits',/remove_all)
;        fits_write, filename, aordark, header ; just use the last header
        
     endfor                     ; for each aorname
  endfor  ; for each exposure time

  ;now what I really want is 3 difference images per position
  ;0.4 - 0.1, 0.1 - 0.02, 0.4 - 0.02
  for p =0, npositions- 1 do begin
     print, 'testing', saveexp(10,10,10,p, 2),saveexp(10,10,10,p, 1)
     diffim = saveexp(*,*,*,p, 2) - saveexp(*,*,*,p,0) ; should be 0.4 - 0.02
     help, diffim
     print, diffim[10,10,5]
     filename = strcompress(dirname+ 'diffdark_'+string(p) +'.fits',/remove_all)
     fits_write, filename, diffim, header  ; just use the last header
  endfor

  end
