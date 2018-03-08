pro make_akaricat
restore, '/Users/jkrick/idlbin/objectnew_akari.sav'

; ok where are they?
;am = where(objectnew.flux11 gt 0 or objectnew.flux15 gt 0 or objectnew.flux18 gt 0)
am = where(objectnew.mips24flux gt 17.3 );and objectnew.acsmag ne 12.828 and objectnew.mips24flux ne 109.63333 and objectnew.mips24flux ne 44.868645)
print, am(1873)
 ;1889 is a bright star that it barfs on when producing the images

openw, outlunred, '/Users/jkrick/akari/allmips.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(am) -1 do  printf, outlunred, 'circle( ', objectnew[am[rc]].ra, objectnew[am[rc]].dec, ' 3")'
close, outlunred
free_lun, outlunred


;print, am
print, 'matched list', n_elements(am)
;print, 'am', am

openw, outlun, '/Users/jkrick/akari/mips24_all_cat.txt', /get_lun
for i = 0, n_elements(am) - 1 do begin

   if objectnew[am(i)].wircjmag gt 0 and objectnew[am(i)].wircjmag ne 99 then begin
      wircjab = 1594.*10^(objectnew[am(i)].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = objectnew[am(i)].wircjmag
   endelse

   if objectnew[am(i)].wirchmag gt 0 and objectnew[am(i)].wirchmag ne 99 then begin
      wirchab = 1024.*10^(objectnew[am(i)].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = objectnew[am(i)].wirchmag
   endelse

   if objectnew[am(i)].wirckmag gt 0 and objectnew[am(i)].wirckmag ne 99 then begin
      wirckab = 666.8*10^(objectnew[am(i)].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = objectnew[am(i)].wirckmag
   endelse

   if objectnew[am(i)].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if objectnew[am(i)].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[am(i)].irac2)
   if objectnew[am(i)].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[am(i)].irac3)
   if objectnew[am(i)].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[am(i)].irac4)
   
   if objectnew[am(i)].imagerr gt 1000. then objectnew[am(i)].imagerr = 1000.
   if objectnew[am(i)].gmagerr gt 1000. then objectnew[am(i)].gmagerr = 1000.
   if objectnew[am(i)].rmagerr gt 1000. then objectnew[am(i)].rmagerr = 1000.
   if objectnew[am(i)].umagerr gt 1000. then objectnew[am(i)].umagerr = 1000.
   
   if objectnew[am(i)].tmassjmag lt 0 then tmassjerr = 99 else tmassjerr = 0.02
   if objectnew[am(i)].tmasshmag lt 0 then tmassherr = 99 else tmassherr = 0.02
   if objectnew[am(i)].tmasskmag lt 0 then tmasskerr = 99 else tmasskerr = 0.02


   printf, outlun, format='(I10,F10.2,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
           am(i), objectnew[am(i)].zphot, objectnew[am(i)].ra,  objectnew[am(i)].dec, objectnew[am(i)].umag,  objectnew[am(i)].umagerr,$
           objectnew[am(i)].gmag, objectnew[am(i)].gmagerr,$
           objectnew[am(i)].rmag, objectnew[am(i)].rmagerr, objectnew[am(i)].acsmag, objectnew[am(i)].acsmagerr,$
           objectnew[am(i)].zmagbest, objectnew[am(i)].zmagerrbest, wircjmagab,  objectnew[am(i)].wircjmagerr, wirchmagab, $
           objectnew[am(i)].wirchmagerr,wirckmagab, objectnew[am(i)].wirckmagerr,$
           objectnew[am(i)].irac1flux,objectnew[am(i)].irac1fluxerr,objectnew[am(i)].irac2flux,objectnew[am(i)].irac2fluxerr,$
           objectnew[am(i)].irac3flux, objectnew[am(i)].irac3fluxerr, objectnew[am(i)].irac4flux, objectnew[am(i)].irac4fluxerr,$
           objectnew[am(i)].flux11, objectnew[am(i)].err11, objectnew[am(i)].flux15, $
           objectnew[am(i)].err15, objectnew[am(i)].flux18, objectnew[am(i)].err18,$
           objectnew[am(i)].mips24flux, objectnew[am(i)].mips24fluxerr,$
           objectnew[am(i)].mips70flux, objectnew[am(i)].mips70fluxerr
 endfor
       
close, outlun
free_lun, outlun
         
plothyperz, am, '/Users/jkrick/akari/mips24_all.ps'


end
