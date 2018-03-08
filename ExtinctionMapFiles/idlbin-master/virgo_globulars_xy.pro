pro virgo_globulars_xy

readcol, '/Users/jkrick/Virgo/IRAC/globulars/ngvs_ra_dec_g24.txt',  ra, dec;, format = '(G10.4, G10.4)'

print, format ='(g10.9)', ra 

fits_read, '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/mosaic_bkgd.fits', data, header

adxy, header, ra, dec, xcenter, ycenter
a = where(xcenter gt 0 and ycenter gt 0 and finite(data(xcenter, ycenter)) eq 1)
xcenter = xcenter(a)
ycenter = ycenter(a)

openw, outlun, '/Users/jkrick/Virgo/IRAC/globulars/ngvs_xy_g24.txt', /get_lun
for r=0, n_elements(xcenter) -1 do  printf, outlun, xcenter(r), '  ', ycenter(r)

close, outlun
free_lun, outlun


end
