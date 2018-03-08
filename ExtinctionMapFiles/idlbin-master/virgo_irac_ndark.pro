pro virgo_irac_ndark

;manually made a median of just the off cluster frames.  
;then read in that median dark
;fits_read, '/Users/jkrick/Virgo/IRAC/r35321856/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
;CD, '/Users/jkrick/Virgo/IRAC/r35321856/ch1/bcd/'      
;command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
;spawn, command1
  
                                ;read in the list of images
;readcol, '/Users/jkrick/Virgo/IRAC/r35321856/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
;for i = 0, n_elements(bcdlist) - 1 do begin
;   fits_read, bcdlist(i), data, header
;   data = data - new_dark
;   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
;   fits_write, outfile, data, header
;endfor
;-------------------------------

fits_read, '/Users/jkrick/Virgo/IRAC/r35322112/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
CD, '/Users/jkrick/Virgo/IRAC/r35322112/ch1/bcd/'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/IRAC/r35322112/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor
;-------------------------------
fits_read, '/Users/jkrick/Virgo/IRAC/r35322368/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
CD, '/Users/jkrick/Virgo/IRAC/r35322368/ch1/bcd/'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/IRAC/r35322368/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor
;-------------------------------
fits_read, '/Users/jkrick/Virgo/IRAC/r35324928/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
CD, '/Users/jkrick/Virgo/IRAC/r35324928/ch1/bcd/'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/IRAC/r35324928/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor
;-------------------------------
fits_read, '/Users/jkrick/Virgo/IRAC/r35325440/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
CD, '/Users/jkrick/Virgo/IRAC/r35325440/ch1/bcd/'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/IRAC/r35325440/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor
;-------------------------------
fits_read, '/Users/jkrick/Virgo/IRAC/r35325952/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
CD, '/Users/jkrick/Virgo/IRAC/r35325952/ch1/bcd/'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/IRAC/r35325952/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor
;-------------------------------
fits_read, '/Users/jkrick/Virgo/IRAC/r35326208/ch1/bcd/med_ch1.fits', new_dark, darkhead

;subtract the dark from all the individual frames.
CD, '/Users/jkrick/Virgo/IRAC/r35326208/ch1/bcd/'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/IRAC/r35326208/ch1/bcd/bcd.txt', bcdlist, format='A', /silent
for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor
;-------------------------------



end
