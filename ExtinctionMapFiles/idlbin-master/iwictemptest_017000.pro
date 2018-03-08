pro iwictemptest

ps_open, filename='/Users/jkrick/iwic/iwic_recovery1/iwictemptest_100_ch1.ps',/portrait,/square,/color
!p.multi=[0,1,2]
;somehow read in fits  files that Patrick sends me
readcol, '/Users/jkrick/iwic/iwic_recovery1/dir_list', dir_names, format="A"

;readcol, '/Users/jkrick/iwic/iwicB/IRAC016400/raw/dir_list', dirnames, format="A",/silent
;readcol, '/Users/jkrick/iracdata/lowrance/IWICA/mipl/dir_list', dirnamesA, format = "A", /silent
;dirnames = [dirnames, dirnamesA]

count = 0
modearr = fltarr(500)
temparr = fltarr(500)
sclkarr = fltarr(500)
hotpixarr = fltarr(500)

biasarr_ch1 = ['B600R3500G3650','B450R3500G3650','B750R3500G3650']
biasarr_ch1 = [biasarr_ch1, biasarr_ch1, biasarr_ch1]
biasarr_ch2 = ['B500R3500G3650','B450R3500G3650','B600R3500G3650']
biasarr_ch2 = [biasarr_ch2, biasarr_ch2, biasarr_ch2]

  for name = 0,   n_elements(dir_names) -1 do begin
     
     CD, dir_names(name)
     spawn, 'pwd'
     command =  ' find . -name "*mipl.fits" >  fits_list.txt'
     spawn, command
     readcol, 'fits_list.txt', fitsnames, format="A", /silent

     for fits = 0, n_elements(fitsnames) - 1 do begin
;        print, 'fitsnames', fitsnames(fits)
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

        
;        if fowlnum eq 8 and waitper eq 44 and ichan eq 2 then begin
        if fowlnum eq 16 and waitper eq 468 and ichan eq 1 then begin
;           print, count, fitsnames(fits)
           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
;           rawdata = rawdata / 12 ; first get it into DN/s
;        rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

;           med = median(rawdata[50:150, 50:150])
           mmm, rawdata, skymod, sigma, skew

     ;then need to write out the fits file
;           newname = strmid(fitsnames(fits), 0,38)
;           newname = strcompress(newname + 'wrap.fits', /remove_all)
;        print, 'would write to ', newname
;           fits_write, newname, rawdata, rawheader

     ;display image on output

;           plotimage, xrange=[1,256],yrange=[1,256], $
;                      bytscl(rawdata, min =min(rawdata), max = max(rawdata) - 4*sigma) ,$ ;skymod+2*sigma
;                      /preserve_aspect, /noaxes,  title= 'raw data'

     ;xyouts the temperature from the headers at that time
;           xyouts, 100 , 30, 'channel '+ string( ichan), charthick = 5
;           xyouts, 100, 50, 'aorID ' + string(aorid), charthick = 5
;           xyouts, 100, 70, 'sclk ' + string(sclk), charthick = 5
;           xyouts, 100, 90,   'AFPAT1B ' + string(AFPAT1B), charthick = 5
;           xyouts, 100, 110,  'AFPAT2B ' + string(AFPAT2B), charthick = 5
;           xyouts, 100, 130,  'AFPAT1E ' + string(AFPAT1E), charthick = 5
;           xyouts, 100, 150,   'AFPAT2E ' + string(AFPAT2E), charthick = 5
;           xyouts, 100, 170, 'mode' + string(skymod), charthick = 5

           ; get the number of hot pixels

           hotpix = where(rawdata gt skymod+(5*sigma))
;           print, 'mode, sigma, hotpix', skymod, sigma, stddev(rawdata), n_elements(hotpix)
;           print, rawdata[hotpix]


           ;save some values for later plotting
           sclkarr(count) = sclk
           modearr(count) = skymod
           temparr(count) = AFPAT2B
           hotpixarr(count) = n_elements(hotpix)
           count = count + 1
           print, 'count', count
        endif

     endfor

  endfor

sclkarr = sclkarr[0:count - 1]
temparr = temparr[0:count - 1]
modearr = modearr[0:count - 1]
hotpixarr = hotpixarr[0:count-1]


plot, temparr, modearr, psym = 2, xtitle = 'Temperature (AFPAT2B)', ytitle = 'Mode DN/S (with sigma clip)', title = '100s  ch 1';, yrange=[300,400], ystyle = 1;, xrange=[44,54], yrange=[0,800]
xyouts, temparr(0)+0.1, modearr(0), biasarr_ch1(0)
xyouts, temparr(5)+0.1, modearr(5), biasarr_ch1(1)
xyouts, temparr(10)+0.1, modearr(10), biasarr_ch1(2)
xyouts, temparr(15)+0.1, modearr(15), biasarr_ch1(0)
xyouts, temparr(20)+0.1, modearr(20), biasarr_ch1(1)
xyouts, temparr(25)+0.1, modearr(25), biasarr_ch1(2)
xyouts, temparr(30)+0.1, modearr(30), biasarr_ch1(0)
xyouts, temparr(35)+0.1, modearr(35), biasarr_ch1(1)
xyouts, temparr(40)+0.1, modearr(40), biasarr_ch1(2)
plot, sclkarr,modearr, psym = 2, xtitle = 'Time (sclk)', ytitle = 'Mode DN/S (with sigma clip)';, yrange=[300,400], ystyle =1;, xrange=[9.2719E8,9.2725E8], yrange=[0,800]
xyouts, sclkarr(0)+0.1, modearr(0), biasarr_ch1(0)
xyouts, sclkarr(5)+0.1, modearr(5), biasarr_ch1(1)
xyouts, sclkarr(10)+0.1, modearr(10), biasarr_ch1(2)
xyouts, sclkarr(15)+0.1, modearr(15), biasarr_ch1(0)
xyouts, sclkarr(20)+0.1, modearr(20), biasarr_ch1(1)
xyouts, sclkarr(25)+0.1, modearr(25), biasarr_ch1(2)
xyouts, sclkarr(30)+0.1, modearr(30), biasarr_ch1(0)
xyouts, sclkarr(35)+0.1, modearr(35), biasarr_ch1(1)
xyouts, sclkarr(40)+0.1, modearr(40), biasarr_ch1(2)
;xyouts, sclkarr(0)+0.1, modearr(0), biasarr_ch1(0)
;xyouts, sclkarr(10)+0.1, modearr(10), biasarr_ch1(1)
;xyouts, sclkarr(20)+0.1, modearr(20), biasarr_ch1(2)
;xyouts, sclkarr(30)+0.1, modearr(30), biasarr_ch1(0)
;xyouts, sclkarr(40)+0.1, modearr(40), biasarr_ch1(1)
;xyouts, sclkarr(50)+0.1, modearr(50), biasarr_ch1(2)
;xyouts, sclkarr(60)+0.1, modearr(60), biasarr_ch1(0)
;xyouts, sclkarr(70)+0.1, modearr(70), biasarr_ch1(1)
;xyouts, sclkarr(80)+0.1, modearr(80), biasarr_ch1(2)

;plot, temparr, hotpixarr, psym = 2, xtitle = 'Temperature (AFPAT2B)', ytitle = 'Number 5sigma pixels';, yrange=[1000,2500], xrange=[44,54]

;print, sorttemp(b)
ps_close, /noprint,/noid

end
