pro subarray

numoobjects = 23*2*47 ; 23 image per AOR, 2 channels, 47 AORs in april calstar pixnoise directory
calstar = replicate({calob, ra:0D,dec:0D,bcdxcen:0D,bcdycen:0D,bcdflux:0D,ch:0I,filename:' ', sclktime:0D},numoobjects)


testname =  '/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/r39167488/ch2/bcd/SPITZER_I2_39167488_0023_0000_2_bcd.fits'
fits_read, testname, data, header
;print, data[20,20,*]

ra =269.84614
dec = 66.049296
get_centroids,testname, t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], RA=ra, DEC=dec;,/silent

for i = 0, n_elements(bcdxcen) - 1 do print, bcdxcen(i), bcdycen(i), bcdflux(i)

end
