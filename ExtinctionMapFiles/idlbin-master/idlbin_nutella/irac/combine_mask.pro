pro combine_mask
;how many background pixels will be left if I combine mask files from
;3 different coordinates?

;central positions of the 3 aors for iwic210 in coords of the irac
;dark field image
x = [1351,1507,1699]
y = [775,971,1219]

fits_read, '/Users/jkrick/IRAC/iwic210/dark_mask.fits', objectmaskdata, objectmaskheader

mask1 = objectmaskdata[1100:1600,525:1025]
mask2 = objectmaskdata[1257:1757,721:1221]
mask3 = objectmaskdata[1450:1950,970:1470]

mask = mask1 + mask2 + mask3

fits_write, '/Users/jkrick/IRAC/iwic210/mask3.fits', mask, objectmaskheader

;how many pixels have zero value? aka are background?

a = where(mask eq 0)
print, 'zero', n_elements(a), ' out of ', 500.*500., ' fraction', n_elements(a) / (500.*500.)
end
