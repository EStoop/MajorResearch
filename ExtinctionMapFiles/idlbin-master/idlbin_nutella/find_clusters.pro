pro findclusters

restore, '/Users/jkrick/idlbin/object.sav'
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');
ch1head = headfits('/Users/jkrick/spitzer/irac/ch1/mosaic.fits')

z1 = 0.0


for zrange = 0, 20 , 1 do begin
 
   z1 =  zrange/10.
   z2 =z1 + 0.2
   
   test = where(object.zphot ge z1 and object.zphot le z2  and object.zphot gt 0.05)
   print, z1,z2,n_elements(test) 
   adxy, ch1head, object[test].ra, object[test].dec, x, y
   
   outfilename = strcompress("/Users/jkrick/nep/clusters/candidate" + string(z1) + string(z2) + ".reg")
   openw, outlun, outfilename, /get_lun
   
   for count = 0, n_elements(test) - 1 do printf, outlun, 'circle(',x(count), y(count), ' 10)'
   
   close, outlun
   free_lun, outlun
   
endfor

end

