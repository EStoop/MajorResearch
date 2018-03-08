pro muxbleed


fits_read, '/Users/jkrick/iwic/iwic_recovery12/IRAC018400/raw/0035749120/IRAC.2.0035749120.0005.0000.01.mipl.fits', rawdata_int, rawheader
barrel = fxpar(rawheader, 'A0741D00')
fowlnum = fxpar(rawheader, 'A0614D00')
pedsig = fxpar(rawheader, 'A0742D00')
ichan = fxpar (rawheader, 'CHNLNUM')
            
            
dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata

fits_write,  '/Users/jkrick/iwic/iwic_recovery12/IRAC018400/raw/0035749120/IRAC.2.0035749120.0005.0000.01.dewrap.fits', rawdata, rawheader

fits_read, '/Users/jkrick/iwic/iwic_recovery12/IRAC018400/bcd/0035749120/IRAC.2.0035749120.0005.0000.1.bcd_fp.fits', bcddata, bcdheader

diff = rawdata / bcddata

fits_write,  '/Users/jkrick/iwic/iwic_recovery12/IRAC018400/diff.fits', diff, rawheader

end
