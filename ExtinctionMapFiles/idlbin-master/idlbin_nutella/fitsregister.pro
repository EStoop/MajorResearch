pro fitsregister

close, /all

;take images which already have wcs in their headers, and put them on the same output pixels.

FITS_READ, '/Users/jkrick/palomar/lfc/coadd_u.fits',udata, uheader
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_g.fits',gdata, gheader
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_r.fits',rdata, rheader
FITS_READ, '/Users/jkrick/palomar/lfc/coadd_i.fits',idata, iheader


HASTROM, udata, uheader, gheader,missing=0
HASTROM, rdata, rheader, gheader,missing=0
HASTROM, idata, iheader, gheader,missing=0

fits_write, '/Users/jkrick/palomar/lfc/u.shift.fits',udata, uheader
fits_write, '/Users/jkrick/palomar/lfc/r.shift.fits',rdata, rheader
fits_write, '/Users/jkrick/palomar/lfc/i.shift.fits',idata, iheader


;then add the frames togehter to make a meaningless image except for object detection
all = udata+rdata+gdata+idata
fits_write, '/Users/jkrick/palomar/lfc/all.shift.fits',all, gheader
end
