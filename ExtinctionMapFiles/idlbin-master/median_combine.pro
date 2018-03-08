pro median_combine, dir_name
 ; dir_name = ['/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037384704/','/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037384960/','/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037385216/','/Users/jkrick/iwic/iwicDAW/IRAC022000/bcd/0037385472/']

;read in the flats for later
 ;fits_read, '/Users/jkrick/irac_pc/pc2/IRAC022300/irac_b1_fa_superskyflat_070709.fits', ch1flat, ch1flathead
 ; fits_read, '/Users/jkrick/irac_pc/pc2/IRAC022300/irac_b2_fa_superskyflat_070709.fits', ch2flat, ch2flathead

;dir_name = ['/Users/jkrick/irac_warm/pc9/IRAC023000/bcd/0038162944/']

;dir_name =  ['/Users/jkrick/Virgo/IRAC/r35321856/ch2/bcd/',  '/Users/jkrick/Virgo/IRAC/r35322112/ch2/bcd/',  '/Users/jkrick/Virgo/IRAC/r35322368/ch2/bcd/',  '/Users/jkrick/Virgo/IRAC/r35324928/ch2/bcd/',  '/Users/jkrick/Virgo/IRAC/r35325440/ch2/bcd/',  '/Users/jkrick/Virgo/IRAC/r35325952/ch2/bcd/',  '/Users/jkrick/Virgo/IRAC/r35326208/ch2/bcd/']

  for j = 0, n_elements(dir_name) - 1 do begin
     CD,dir_name(j)

     command1 =  ' ls *_bcd.fits >  bcd_list.txt'
     command2 = 'cat bcd_list.txt | grep SPITZER_I1 >bcd_ch1.txt'
;     command2 = 'cat bcd_list.txt | grep SPITZER_I2 > bcd_ch2.txt'
     
     commands = [command1, command2];, command3]
     for i = 0, n_elements(commands) -1 do spawn, commands(i)
     
     readcol,dir_name(j) + 'bcd_ch1.txt', bcdch1list, format="A", /silent
;     readcol,dir_name(j) + 'bcd_ch2.txt', bcdch1list, format="A", /silent

     print, 'n ch1list,', n_elements(bcdch1list);, n_elements(bcdch2list)
     bigim = fltarr(32, 32, 128, n_elements(bcdch1list));n_elements(bcdch1list))
;     bigim = fltarr(256, 256, 70);n_elements(bcdch1list))
     c = 0
        
     for i =0, n_elements(bcdch1list) -1 do begin
;        print, i, bcdch1list(i), c
;     for i =0, 70 - 1 do begin
;     for i =n_elements(bcdch1list) -70, n_elements(bcdch1list) -1 do begin
        fits_read, bcdch1list(i), im, imheader
        ;im = im * ch1flat
        ;im = im + 0.25
        ;im = im / ch1flat
;        bigim(0,0,c) = im
        bigim(0,0,0,c) = im
        c = c + 1
     endfor
     
;     medarr, bigim, meddark
     meddark = median(bigim, dimension = 4)
     fits_write, dir_name(j) + 'med_ch1.fits', meddark, imheader
;     fits_write, dir_name(j) + 'med_ch1_beg.fits', meddark, imheader
;     fits_write, dir_name(j) + 'med_ch1_end.fits', meddark, imheader
;----

;     bigim = fltarr(256, 256, n_elements(bcdch2list) - 58)
;     for i = 59, n_elements(bcdch2list) - 1 do begin
;        fits_read, bcdch2list(i), im, imheader
;        bigim(0,0,i-59) = im
;        ;im = im * ch2flat
;        ;im = im + 0.25
;        ;im = im / ch2flat
;        ;bigim(0,0,i) = im
       
;     endfor
     
;     medarr, bigim, meddark
     
;     fits_write, dir_name(j) + 'med_ch2.fits', meddark, imheader
  endfor

end
