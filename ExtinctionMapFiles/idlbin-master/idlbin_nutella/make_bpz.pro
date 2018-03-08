pro make_bpz


close, /all

restore, 'object.sav'

openw, outlun, "/Users/jkrick/bin/bpz.1.98b/bpz_cat.txt",/get_lun
count = 1

for num = 0, n_elements(object.uxcenter) -1 do begin

if object[num].irac1 lt 0 then mag1 = -99. else mag1 = 6.12 - 2.5*alog10(1E-6*object[num].irac1) 
if object[num].irac2 lt 0 then mag2 = -99. else mag2 = 5.64 - 2.5*alog10(1E-6*object[num].irac2) 
if object[num].irac3 lt 0 then mag3 = -99. else mag3 = 5.16 - 2.5*alog10(1E-6*object[num].irac3)
if object[num].irac4 lt 0 then mag4 = -99. else mag4 = 4.52 - 2.5*alog10(1E-6*object[num].irac4)
if object[num].irac1 lt 0 then err1 = 26. else err1 =  0.05*(6.12 - 2.5*alog10(1E-6*object[num].irac1))
if object[num].irac2 lt 0 then err2 = 26. else err2 =  0.05*(5.64 - 2.5*alog10(1E-6*object[num].irac2))
if object[num].irac3 lt 0 then err3 = 26. else err3 =  0.05*(5.16 - 2.5*alog10(1E-6*object[num].irac3))
if object[num].irac4 lt 0 then err4 = 26. else err4 =  0.05*(4.52 - 2.5*alog10(1E-6*object[num].irac4))

;reject some stars based on fwhm's
   if object[num].gfwhm gt 0 and object[num].gfwhm lT 5.6  then begin
;      print, object[num].gfwhm
      ;don't want these stars
   endif else begin
      printf, outlun, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',num+1, object[num].umaga, object[num].gmaga, object[num].rmaga, $
              object[num].imaga, mag1,mag2,mag3,mag4,  object[num].umagerra, object[num].gmagerra, $
              object[num].rmagerra, object[num].imagerra, err1,err2,err3,err4

      count = count + 1

   endelse


endfor

close, outlun
free_lun, outlun

end

;reject objects which only have 1 flux measurement
;isn't working since there are some bright objects which have 2 measurements
;reject these based on bad probabilities after running hyperz

;   sum = object[num].umaga+ object[num].gmaga+ object[num].rmaga+ object[num].imaga+ $
;         mag1+mag2+mag3+mag4

;   if sum LT .0 then begin

;      print, count, object[num].umaga, object[num].gmaga, object[num].rmaga, object[num].imaga, $
;             mag1,mag2,mag3,mag4 ,  $
;             object[num].umagerra, object[num].gmagerra, object[num].rmagerra, object[num].imagerra, $
;            err1,err2,err3,err4

;   endif
