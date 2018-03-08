PRO count_flux
close,/all

fits_read, "/n/Godiva4/jkrick/A2734/original/test.fits", data, header
inner = 0.
innersum = 0.D
outer = 0.
outersum = 0.
sum = 0.D
npix = 0.D
npix2 = 0.D
true = 0
FOR x = 2400, 2500,1 DO BEGIN
    FOR y= 1, 100, 1 DO BEGIN
        innersum = innersum + data[x,y]
    ENDFOR
ENDFOR

print,innersum
END
