pro cobe_dirbe
; ftab_help, '/Users/jkrick/irac_darks/DIRBE/DIRBE_CIO_P3B_90159.FITS'
;figure out which DIRBE pixel to use
;  pixnum = coorconv([270.,66.56], infmt = 'H',  outfmt = 'P', inco = 'Q', outco = 'R9') ;the NEP
;  pixnum = coorconv([18.,66.56], infmt = 'H',  outfmt = 'P', inco = 'Q', outco = 'R9') ;the NEP
pixnum = coorconv([17.6667,69.], infmt = 'H',  outfmt = 'P', inco = 'Q', outco = 'R9') ; the dark field center
pixnum = pixnum[0]
;pixnum = 23454
help, pixnum

  dirbedir =  '/Users/jkrick/irac_darks/DIRBE/'
  indexdir =  '/Users/jkrick/irac_darks/DIRBE/indx'
  
  CD, dirbedir                   
  command  =  " ls *.FITS > /Users/jkrick/irac_darks/DIRBE/indx/indexlist.txt"  
  spawn, command
  
  readcol,'/Users/jkrick/irac_darks/DIRBE/indx/indexlist.txt',indexname, format = 'A', /silent
  timearr = fltarr(285*10)
  phot3Aarr = fltarr(285*10)
  phot3barr = fltarr(285*10)
  phot3carr = fltarr(285*10)
  phot4arr = fltarr(285*10)
  phot5arr = fltarr(285*10)
  testarr =  fltarr(285*10)
  moon2losarr =  fltarr(285*10)
  jup2losarr =  fltarr(285*10)
  xsnoisearr =  fltarr(285*10)
  c = 0
  for i = 0, n_elements(indexname) - 1 do begin
     ftab_ext, indexname(i), 'pixel_no,phot3A,Time,phot3B,phot3C,phot04,phot05,moon2los,jup2los,xsnoise', pixel_no,phot3A,Time,phot3b, phot3c, phot04, phot05,moon2los,jup2los,xsnoise
     a = where(pixel_no eq pixnum, count)
     if count gt 0 then begin
        print,indexname(i), count
        for t = 0, count - 1 do begin
           print, pixnum, pixel_no[a[t]]
           timearr[c] = Time[a[t]]
           phot3aarr[c] = phot3a[a[t]]
            phot3barr[c] = phot3b[a[t]]
           phot3carr[c] = phot3c[a[t]]
           phot4arr[c] = phot04[a[t]]
           phot5arr[c] = phot05[a[t]]
          testarr[c] = pixel_no[a[t]]
          moon2losarr[c] = moon2los[a[t]]
          jup2losarr[c] = jup2los[a[t]]
           xsnoisearr[c] = xsnoise[a[t]]
          c = c + 1
        endfor

     endif

  endfor

timearr = timearr[0:c-1]
phot3aarr = phot3aarr[0:c-1]
phot3barr = phot3barr[0:c-1]
phot3carr = phot3carr[0:c-1]
phot4arr = phot4arr[0:c-1]
phot5arr = phot5arr[0:c-1]
testarr = testarr[0:c-1]
moon2losarr = moon2losarr[0:c-1]
jup2losarr = jup2losarr[0:c-1]
xsnoisearr = xsnoisearr[0:c-1]

print, 'test', mean(testarr), stddev(testarr)
print, 'time', timearr[0:10]
;p = plot(timearr, phot3aarr, '1*')

  save, /all, filename='/Users/jkrick/irac_darks/dirbe/dirbe.sav'

end
