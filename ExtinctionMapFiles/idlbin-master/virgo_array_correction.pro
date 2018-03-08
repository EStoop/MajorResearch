pro virgo_array_correction
!P.multi = [0,2,1]

  dirname = '/Users/jkrick/Virgo/IRAC/'
;------------------------------------------------
;ch1
;read in the array location depenedence images for the warm mission
  fits_read, '/Users/jkrick/virgo/irac/array_location/ch1_photcorr_ap_5.fits', ch1arraydata, ch1arrayheader

  ;read in the list of all the cbcd's used in the mosaic
  readcol, '/Users/jkrick/virgo/irac/cbcd_ch1_list_2.txt', filename,format="A"
  filename = strcompress(dirname + strmid(filename, 2),/remove_all)
  ;print, filename

  ;make a new list of filenames to be used for the output.
  newfilename = strcompress(strmid(filename, 0,79) + '_array.fits',/remove_all)
  ;print, newfilename

  ;make an output list of newfilenames as input to imcombine
  openw, outlunred, '/Users/jkrick/virgo/irac/cbcd_array_ch1_list.txt', /get_lun
  for rc=0, n_elements(filename) -1 do  printf, outlunred, newfilename(rc)
  close, outlunred
  free_lun, outlunred
  
;for i = 0, n_elements(filename) - 1 do begin
;   fits_read, filename(i), data, header
;   fits_write, newfilename(i), ch1arraydata, header
;endfor
  
;now read in the combined array location dependence image and do some statistics
fits_read, '/Users/jkrick/Virgo/IRAC/array_location/ch1/ch1_Combine-mosaic/mosaic.fits', data, header
plothist, data, xarr, yarr, bin = 0.01, xrange = [0.9, 1.0],/ylog, xtitle = 'Array dependent correction (%)', ytitle = 'Number', title = 'CH1'


;------------------------------------------------
;ch2
;read in the array location depenedence images for the warm mission
  fits_read, '/Users/jkrick/virgo/irac/array_location/ch2_photcorr_ap_5.fits', ch2arraydata, ch2arrayheader

  ;read in the list of all the cbcd's used in the mosaic
  readcol, '/Users/jkrick/virgo/irac/cbcd_ch2_list_2.txt', filename,format="A"
  filename = strcompress(dirname + strmid(filename, 2),/remove_all)
  ;print, filename

  ;make a new list of filenames to be used for the output.
  newfilename = strcompress(strmid(filename, 0,79) + '_array.fits',/remove_all)
  ;print, newfilename

  ;make an output list of newfilenames as input to imcombine
  openw, outlunred, '/Users/jkrick/virgo/irac/cbcd_array_ch2_list.txt', /get_lun
  for rc=0, n_elements(filename) -1 do  printf, outlunred, newfilename(rc)
  close, outlunred
  free_lun, outlunred
  
;for i = 0, n_elements(filename) - 1 do begin
;   fits_read, filename(i), data, header
;   fits_write, newfilename(i), ch2arraydata, header
;endfor
  
;now read in the combined array location dependence image and do some statistics
fits_read, '/Users/jkrick/Virgo/IRAC/array_location/ch2/ch2_Combine-mosaic/mosaic.fits', data, header
plothist, data, xarr, yarr, bin = 0.01, xrange = [0.9, 1.0],/ylog, xtitle = 'Array dependent correction (%)', ytitle = 'Number', title = 'CH2'

end
