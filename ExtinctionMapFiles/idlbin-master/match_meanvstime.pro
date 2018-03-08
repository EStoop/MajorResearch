pro match_meanvstime
;taking the information from fit_meanvstime.pro
;make a correction to both the ramp and the background levels within each tile.
restore,'/Users/jkrick/virgo/irac/fitmean_ch2.sav'
ch = 1
  for a =0 ,n_elements(sortaor) - 1 do begin
   
      cd, strcompress(dirloc + sortaor(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
      spawn, 'pwd'
      ;command =  "find . -name 'corr_SPITZER_I1*0_ndark_ff.fits' > ndark_ff_ch1.txt"
      command =  "find . -name 'corr_SPITZER_I2*0_ndark.fits' > ndark_ch2.txt"
      spawn, command
      
      ;read in a list of all ndark and ff corrected images
      readcol,strcompress( 'ndark_ff_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
      ;readcol,strcompress( 'ndark_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
      ;print, 'would be adding', 0.28 - levelarr(a)
      print, 'would be adding', 0.16 - levelarr(a)

      for i =1, n_elements(fitsname) - 1 do begin
         fits_read, fitsname(i), data, header
         ;first work on the overall levels
;         data = data +( 0.28 - levelarr(a))
         ; work on the ramp and overall level
        ;data = data + ( 0.28 - (levelarr(a) - result1arr(a)*exp(-result2arr(a)*i)))
         data = data + ( 0.16 - (levelarr(a) - result1arr(a)*exp(-result2arr(a)*i)))

        newname = strcompress(strmid(fitsname(i), 0, 37) + 'ndark_ff_match.fits',/remove_all)
;        newname = strcompress(strmid(fitsname(i), 0, 37) + 'ndark_match.fits',/remove_all)
         fits_write, newname, data, header
      endfor  ;for each fits file
   endfor ; for each AOR




end
