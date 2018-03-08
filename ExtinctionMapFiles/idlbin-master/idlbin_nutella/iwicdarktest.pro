pro iwicdarktest
;measure the true dark current in our 50K IWICA & B data from first IWIC
;still need to add IWIC B AORs and 6, and 200s images
;still need to change the last part to measure dark between 100 and 200s
;not between 12 and 100
;the expected value at 50K is ~500-600 e/s
;caveats are saturation and therefore a different gain


;ps_open, filename='/Users/jkrick/iwic/iwicdarktest_ch1.ps',/portrait,/square,/color
!p.multi=[0,1,1]
;somehow read in fits  files that Patrick sends me
readcol, '/Users/jkrick/iwic/iwicB/IRAC016400/raw/dir_list', dirnames, format="A",/silent
readcol, '/Users/jkrick/iwic/iwicA/IRAC016300/raw/dir_list_short', dirnamesA, format = "A", /silent
;dirnames = [dirnames, dirnamesA]
dirnames = dirnamesA ;not using B for the moment want to just use setpoint tests in B
gain = [3.3,3.3,3.3,3.3,3.3,3.3]  ; change this for adding B

count2 = 0
modearr_2 = fltarr(500)
temparr_2 = fltarr(500)
sclkarr_2 = fltarr(500)
hotpixarr_2 = fltarr(500)

count6 = 0
modearr_6 = fltarr(500)
temparr_6 = fltarr(500)
sclkarr_6 = fltarr(500)
hotpixarr_6 = fltarr(500)

count12 = 0
modearr_12 = fltarr(500)
temparr_12 = fltarr(500)
sclkarr_12 = fltarr(500)
hotpixarr_12 = fltarr(500)

count100 = 0
modearr_100 = fltarr(500)
temparr_100 = fltarr(500)
sclkarr_100 = fltarr(500)
hotpixarr_100 = fltarr(500)

count200 = 0
modearr_200 = fltarr(500)
temparr_200 = fltarr(500)
sclkarr_200 = fltarr(500)
hotpixarr_200 = fltarr(500)

  for name = 0,   n_elements(dirnames) -1 do begin
     
     CD, dirnames(name)

     command =  ' find . -name "*mipl.fits" >  fits_list.txt'
     spawn, command
     readcol, 'fits_list.txt', fitsnames, format="A", /silent

     for fits = 0, n_elements(fitsnames) - 1 do begin
        fits_read, fitsnames(fits), rawdata_int, rawheader

        barrel = fxpar(rawheader, 'A0741D00')
        fowlnum = fxpar(rawheader, 'A0614D00')
        waitper = fxpar(rawheader, 'A0615D00')
        pedsig = fxpar(rawheader, 'A0742D00')
        ichan = fxpar (rawheader, 'CHNLNUM')
        AFPAT1B = fxpar(rawheader, 'A0653E00')
        AFPAT2B = fxpar(rawheader, 'A0655E00')
        AFPAT1E = fxpar(rawheader, 'A0662E00')
        AFPAT2E = fxpar(rawheader, 'A0664E00')
        sclk = fxpar(rawheader, 'A0601D00')
        aorid = fxpar(rawheader, 'A0553D00')

        
        if fowlnum eq 4 and waitper eq 2 and ichan eq 1 then begin  ;2s
;           print, count, fitsnames(fits)

           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
;           rawdata = rawdata / 2 ; first get it into DN/s
;        rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

;           med = median(rawdata[50:150, 50:150])
           mmm, rawdata, skymod, sigma, skew

           ;save some values for later plotting
           sclkarr_2(count2) = sclk
           modearr_2(count2) = skymod * gain(name) ; save in electrons
           temparr_2(count2) = AFPAT2B
           hotpixarr_2(count2) = n_elements(hotpix)
           count2 = count2 + 1

        endif
        if fowlnum eq 8 and waitper eq 14 and ichan eq 1 then begin ;6s
;           print, count, fitsnames(fits)

           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
;           rawdata = rawdata / 12 ; first get it into DN/s
;        rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

;           med = median(rawdata[50:150, 50:150])
           mmm, rawdata, skymod, sigma, skew

           ;save some values for later plotting
           sclkarr_6(count6) = sclk
           modearr_6(count6) = skymod * gain(name) ; save in electrons
           temparr_6(count6) = AFPAT2B
           hotpixarr_6(count6) = n_elements(hotpix)
           count6 = count6 + 1

        endif

        if fowlnum eq 8 and waitper eq 44 and ichan eq 1 then begin  ;12s
;           print, count, fitsnames(fits)

           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
;           rawdata = rawdata / 12 ; first get it into DN/s
;        rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

;           med = median(rawdata[50:150, 50:150])
           mmm, rawdata, skymod, sigma, skew

           ;save some values for later plotting
           sclkarr_12(count12) = sclk
           modearr_12(count12) = skymod * gain(name) ; save in electrons
           temparr_12(count12) = AFPAT2B
           hotpixarr_12(count12) = n_elements(hotpix)
           count12 = count12 + 1

        endif
        if fowlnum eq 16 and waitper eq 468 and ichan eq 1 then begin ;100s
;           print, count, fitsnames(fits)

           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
;           rawdata = rawdata / 12 ; first get it into DN/s
;        rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr
           newname = strmid(fitsnames(fits), 0,38)
           newname = strcompress(newname + 'wrap.fits', /remove_all)
;           print, 'would write to ', newname
           fits_write, newname, rawdata, rawheader

;           med = median(rawdata[50:150, 50:150])
           mmm, rawdata, skymod, sigma, skew

           ;save some values for later plotting
           sclkarr_100(count100) = sclk
           modearr_100(count100) = skymod * gain(name) ; save in electrons
           temparr_100(count100) = AFPAT2B
           hotpixarr_100(count100) = n_elements(hotpix)
           count100 = count100 + 1

        endif
        if fowlnum eq 32  and waitper eq 936 and ichan eq 1 then begin ;200s
           print, count200, fitsnames(fits)
;           spawn, 'pwd'
           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
;           rawdata = rawdata / 12 ; first get it into DN/s
;        rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr
           newname = strmid(fitsnames(fits), 0,38)
           newname = strcompress(newname + 'wrap.fits', /remove_all)
;           print, 'would write to ', newname
           fits_write, newname, rawdata, rawheader

;           med = median(rawdata[50:150, 50:150])
           mmm, rawdata, skymod, sigma, skew

           ;save some values for later plotting
           sclkarr_200(count200) = sclk
           modearr_200(count200) = skymod * gain(name) ; save in electrons
           temparr_200(count200) = AFPAT2B
           hotpixarr_200(count200) = n_elements(hotpix)
           count200 = count200 + 1

        endif

     endfor ;for all fits files within that directory

  endfor  ;for all directories for iwicA and iwicB           
print, 'count2', count2
print, 'cout 6', count6
print, 'count 12', count12
print, 'count 100', count100
print, 'count200', count200

  sclkarr_2 = sclkarr_2[0:count2 - 1]
  temparr_2 = temparr_2[0:count2 - 1]
  modearr_2 = modearr_2[0:count2 - 1]
  hotpixarr_2 = hotpixarr_2[0:count2-1]

  sclkarr_6 = sclkarr_6[0:count6 - 1]
  temparr_6 = temparr_6[0:count6 - 1]
  modearr_6 = modearr_6[0:count6 - 1]
  hotpixarr_6 = hotpixarr_6[0:count6-1]

  sclkarr_12 = sclkarr_12[0:count12 - 1]
  temparr_12 = temparr_12[0:count12 - 1]
  modearr_12 = modearr_12[0:count12 - 1]
  hotpixarr_12 = hotpixarr_12[0:count12-1]

  sclkarr_100 = sclkarr_100[0:count100 - 1]
  temparr_100 = temparr_100[0:count100 - 1]
  modearr_100 = modearr_100[0:count100 - 1]
  hotpixarr_100 = hotpixarr_100[0:count100-1]

  sclkarr_200 = sclkarr_200[0:count200 - 1]
  temparr_200 = temparr_200[0:count200 - 1]
  modearr_200 = modearr_200[0:count200 - 1]
  hotpixarr_200 = hotpixarr_200[0:count200-1]

;want an overall mean(mode) value for DN at each exptime
  average_12 = mean(modearr_12)
  average_2 = mean(modearr_2)
  average_100 = mean(modearr_100)
  average_6 = mean(modearr_6)
  average_200 = mean(modearr_200)

  x = [2,6,12,100,200]
  y = [average_2, average_6,average_12, average_100, average_200]
  plot, x, y, psym = 2, title = '~50K', xrange=[0,200], ytitle = 'electrons', xtitle = 'exposure time'


  ;assume we can make a measurement between 2 and 12
  ;XXX will want to change this when we get all exptimes

  delta_e = average_12 - average_6
  delta_s = 12 - 6

  dark = delta_e / delta_s

  ;subtract off the background value we know in the dark field
  dark = dark - 2.  ;2 e/s in ch1
  print, 'dark in e/s at 50K', dark

;ps_close, /noprint,/noid

end
