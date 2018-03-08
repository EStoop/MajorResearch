pro A1914



CD, '/Users/jkrick/Virgo/A1914/r31483136/ch1/bcd'      
command1 =  ' ls *_cbcd*.fits >  bcd.txt'   ;list all images in that directory
spawn, command1
  
                                ;read in the list of images
readcol, '/Users/jkrick/Virgo/A1914/r31483136/ch1/bcd/bcd.txt', bcdlist, format='A', /silent

;median combine everything to make my own dark image
;median_combine,  '/Users/jkrick/Virgo/A1914/r31483136/ch1/bcd/'

;actually manyally made a median of just the off cluster frames.  

fits_read, '/Users/jkrick/Virgo/A1914/r31483136/ch1/bcd/med_ch1_off.fits', new_dark, darkhead

;subtract the dark from all the individual frames.

for i = 0, n_elements(bcdlist) - 1 do begin
   fits_read, bcdlist(i), data, header
   data = data - new_dark
   outfile = strcompress(strmid(bcdlist(I), 0,31)+ '_ndark.fits', /remove_all)
   fits_write, outfile, data, header
endfor


;what is the background subtraction doing.
;fits_read, '/Users/jkrick/virgo/A1914/r31483136/ch1/bcd/SPITZER_I1_31483136_0167_0000_1_ndark.fits', ndark, ndheader
;fits_read, '/Users/jkrick/virgo/A1914/r31483136/ch1/bcd/Medfilter-mosaic/SPITZER_I1_31483136_0167_0000_1_ndark_minback.fits', mindark, minheader

;back = ndark - mindark
;fits_write, '/Users/jkrick/virgo/A1914/r31483136/ch1/bcd/back_0167-11.fits', back, minheader

end
