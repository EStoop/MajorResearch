pro mips24_ra


fits_read, '/Users/jkrick/spitzer/mips/mips24/dao/newastrometry/mosaic.fits', data, header

readlargecol, "/Users/jkrick/Spitzer/mips/mips24/dao/mips24.phot", num,  xcen,ycen, mag, magerr, back, niter, sh, junk1, junk2, format = 'A'

xyad, header, xcen, ycen, ra, dec

openw, outlun, '/Users/jkrick/spitzer/mips/mips24/dao/newastrometry/mips24.ra.phot', /get_lun
for i = 0, n_elements(xcen) - 1 do printf, outlun, format='(I10,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', num(i),ra(i) ,dec(i), xcen(i),ycen(i), mag(i), magerr(i), back(i), niter(i), sh(i), junk1(i), junk2(i)

close, outlun
free_lun, outlun

end
