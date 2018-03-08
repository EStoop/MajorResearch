pro display_raw_full

fits_read, '/Users/jkrick/irac_warm/hatp7/r28311040/ch2/raw/SPITZER_I2_0028311040_0007_0001_01_dce.fits', data, rawheader
;mc = reform(data, 32, 32, 64)
;tvscl, mc(*,*,0)

;img = data
; convert to real
;img=img*1.

; fix the problem with unsigned ints
;fix=where((img LT 0),count)
;if (count GT 0) then img[fix]=img[fix]+65536.

; flip the InSb
;img=65535.-img

;fits_write, '/Users/jkrick/irac_warm/snapshots/r44497152/ch2/raw/test_0023.fits', img, header
;then use ds9 to display

barrel = fxpar(rawheader, 'A0741D00')
fowlnum = fxpar(rawheader, 'A0614D00')
pedsig = fxpar(rawheader, 'A0742D00')
ichan = fxpar (rawheader, 'CHNLNUM')

dewrap2, data, ichan, barrel, fowlnum, pedsig, 0, rawdata

print, rawdata[3:6,87:91]
print, 'max', max(rawdata[1:20, 50:150])
end
