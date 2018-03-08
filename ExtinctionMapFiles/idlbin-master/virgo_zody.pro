pro virgo_zody

readcol, '/Users/jkrick/Virgo/IRAC/s18p14/ch1_zody.txt', filename, zody_est,format='(A,F10.8)'
print, 'ch1', min(zody_est), max(zody_est)

readcol, '/Users/jkrick/Virgo/IRAC/s18p14/ch2_zody.txt', filename, zody_est,format='(A,F10.8)'
print, 'ch2', min(zody_est), max(zody_est)


end
