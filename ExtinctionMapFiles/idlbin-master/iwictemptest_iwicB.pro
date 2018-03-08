pro iwictemptest

ps_open, filename='/Users/jkrick/iwic/iwictemptest_ch1_all.ps',/portrait,/square,/color
!p.multi=[0,1,2]
;somehow read in fits  files that Patrick sends me
readcol, '/Users/jkrick/iwic/iwicB/IRAC016400/raw/dir_list', dirnames, format="A",/silent
readcol, '/Users/jkrick/iracdata/lowrance/IWICA/mipl/dir_list', dirnamesA, format = "A", /silent
dirnames = [dirnames, dirnamesA]

count = 0
modearr = fltarr(500)
temparr = fltarr(500)
sclkarr = fltarr(500)
hotpixarr = fltarr(500)

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

        
        if fowlnum eq 4 and waitper eq 2 and ichan eq 1 then begin
;           print, count, fitsnames(fits)

           dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
           rawdata = rawdata / 12 ; first get it into DN/s
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

        endif

     endfor

  endfor

sclkarr = sclkarr[0:count - 1]
temparr = temparr[0:count - 1]
modearr = modearr[0:count - 1]
hotpixarr = hotpixarr[0:count-1]

;want different symbols for IWICA and IWICB
sortsclk = sclkarr[sort(sclkarr)]
sorttemp = temparr[sort(sclkarr)]
sortmode = modearr[sort(sclkarr)]
sorthot = hotpixarr[sort(sclkarr)]

a = where(sortsclk lt 9.2723E8)
b = where(sortsclk gt 9.2723E8)

plot, sorttemp(a), sortmode(a), psym = 2, xtitle = 'Temperature (AFPAT2B)', ytitle = 'Mode DN/S (with sigma clip)', xrange=[44,54], yrange=[0,800]
oplot, sorttemp(b), sortmode(b), psym = 6

plot, sortsclk(a), sortmode(a), psym = 2, xtitle = 'Time (sclk)', ytitle = 'Mode DN/S (with sigma clip)', xrange=[9.2719E8,9.2725E8], yrange=[0,800]
oplot, sortsclk(b), sortmode(b), psym = 6

plot, sorttemp(a), sorthot(a), psym = 2, xtitle = 'Temperature (AFPAT2B)', ytitle = 'Number 5sigma pixels', yrange=[1000,2500], xrange=[44,54]
oplot, sorttemp(b), sorthot(b), psym = 6

print, sorttemp(b)
ps_close, /noprint,/noid

end
