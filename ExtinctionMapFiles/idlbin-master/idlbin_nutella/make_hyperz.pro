pro make_hyperz


close, /all

restore, '/Users/jkrick/idlbin/objectnew.sav'


openw, outlun, "/Users/jkrick/ZPHOT/hyperz_cat.txt",/get_lun
count = 1l
starcount = 0
for num = 0l,  n_elements(objectnew.ura) -1 do begin
;AB magnitudes
   ;don't convert things which are already -1, or -99
   
   if objectnew[num].flamjmag gt 0 and objectnew[num].flamjmag ne 99 then begin
      fab = 1594.*10^(objectnew[num].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = objectnew[num].flamjmag
   endelse

   if objectnew[num].wircjmag gt 0 and objectnew[num].wircjmag ne 99 then begin
      wircjab = 1594.*10^(objectnew[num].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = objectnew[num].wircjmag
   endelse

   if objectnew[num].wirchmag gt 0 and objectnew[num].wirchmag ne 99 then begin
      wirchab = 1024.*10^(objectnew[num].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = objectnew[num].wirchmag
   endelse

   if objectnew[num].wirckmag gt 0 and objectnew[num].wirckmag ne 99 then begin
      wirckab = 666.8*10^(objectnew[num].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = objectnew[num].wirckmag
   endelse

   if objectnew[num].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if objectnew[num].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[num].irac2)
   if objectnew[num].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[num].irac3)
   if objectnew[num].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*objectnew[num].irac4)
   
   if objectnew[num].imagerr gt 1000. then objectnew[num].imagerr = 1000.
   if objectnew[num].gmagerr gt 1000. then objectnew[num].gmagerr = 1000.
   if objectnew[num].rmagerr gt 1000. then objectnew[num].rmagerr = 1000.
   if objectnew[num].umagerr gt 1000. then objectnew[num].umagerr = 1000.
   if objectnew[num].zmagerrbest gt 1000. then objectnew[num].zmagerrbest = 1000.
   
   if objectnew[num].tmassjmag lt 0 then tmassjerr = 99 else tmassjerr = 0.02
   if objectnew[num].tmasshmag lt 0 then tmassherr = 99 else tmassherr = 0.02
   if objectnew[num].tmasskmag lt 0 then tmasskerr = 99 else tmasskerr = 0.02

   nflux = 0
   if objectnew[num].umag gt 0 and objectnew[num].umag ne 99 then nflux = nflux + 1
   if objectnew[num].gmag gt 0 and objectnew[num].gmag ne 99 then nflux = nflux + 1
   if objectnew[num].rmag gt 0 and objectnew[num].rmag ne 99 then nflux = nflux + 1
   if objectnew[num].imag gt 0 and objectnew[num].imag ne 99 then nflux = nflux + 1
   if objectnew[num].acsmag gt 0 and objectnew[num].acsmag ne 99 then nflux=nflux+1
   if objectnew[num].zmagbest gt 0 and objectnew[num].zmagbest ne 99 then nflux = nflux +1
   if objectnew[num].flamjmag gt 0 and objectnew[num].flamjmag ne 99 then nflux = nflux + 1
   if objectnew[num].wircjmag gt 0 and objectnew[num].wircjmag ne 99 then nflux = nflux + 1
   if objectnew[num].wirchmag gt 0 and objectnew[num].wirchmag ne 99 then nflux = nflux + 1
   if objectnew[num].wirckmag gt 0 and objectnew[num].wirckmag ne 99 then nflux = nflux + 1
   if objectnew[num].tmassjmag gt 0 and objectnew[num].tmassjmag ne 99 then nflux = nflux + 1
   if objectnew[num].tmasshmag gt 0 and objectnew[num].tmasshmag ne 99 then nflux = nflux + 1
   if objectnew[num].tmasskmag gt 0 and objectnew[num].tmasskmag ne 99 then nflux = nflux + 1
   if objectnew[num].irac1mag gt 0 and objectnew[num].irac1mag ne 99 then nflux = nflux + 1
   if objectnew[num].irac2mag gt 0 and objectnew[num].irac2mag ne 99 then nflux = nflux + 1
   if objectnew[num].irac3mag gt 0 and objectnew[num].irac3mag ne 99 then nflux = nflux + 1
   if objectnew[num].irac4mag gt 0 and objectnew[num].irac4mag ne 99 then nflux = nflux + 1
   if objectnew[num].mips24mag gt 0 and objectnew[num].mips24mag ne 99 then nflux = nflux + 1

   if nflux gt 0 then  printf, outlun, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
                 num, objectnew[num].umag , objectnew[num].gmag, objectnew[num].rmag, $
                 objectnew[num].imag,  objectnew[num].acsmag, objectnew[num].zmagbest, jmagab, wircjmagab, $
                 objectnew[num].tmassjmag,  wirchmagab, objectnew[num].tmasshmag, wirckmagab, $
                 objectnew[num].tmasskmag, objectnew[num].irac1mag,objectnew[num].irac2mag,$
                 objectnew[num].irac3mag,objectnew[num].irac4mag,   $
                 objectnew[num].umagerr, objectnew[num].gmagerr, $
                 objectnew[num].rmagerr, objectnew[num].imagerr, objectnew[num].acsmagerr, objectnew[num].zmagerrbest, objectnew[num].flamjmagerr, $
                 objectnew[num].wircjmagerr,tmassjerr, objectnew[num].wirchmagerr, tmassherr, objectnew[num].wirckmagerr,$
                 tmasskerr,objectnew[num].irac1magerr,objectnew[num].irac2magerr,$
                 objectnew[num].irac3magerr, objectnew[num].irac4magerr;err1,err2,err3,err4
;took out acs until photometry gets sorted out.

         count = count + 1
  

endfor
print, "count", count
print, "num", num
close, outlun
free_lun, outlun

end

;reject objectnews which only have 1 flux measurement
;isn't working since there are some bright objectnews which have 2 measurements
;reject these based on bad probabilities after running hyperz

;   sum = objectnew[num].umaga+ objectnew[num].gmaga+ objectnew[num].rmaga+ objectnew[num].imaga+ $
;         mag1+mag2+mag3+mag4

;   if sum LT .0 then begin

;      print, count, objectnew[num].umaga, objectnew[num].gmaga, objectnew[num].rmaga, objectnew[num].imaga, $
;             mag1,mag2,mag3,mag4 ,  $
;             objectnew[num].umagerra, objectnew[num].gmagerra, objectnew[num].rmagerra, objectnew[num].imagerra, $
;            err1,err2,err3,err4

;   endif

;vega magnitudes

;   if objectnew[num].irac1 lt 0 then mag1 = -1. else mag1 = 6.12 - 2.5*alog10(1E-6*objectnew[num].irac1) 
;   if objectnew[num].irac2 lt 0 then mag2 = -1. else mag2 = 5.64 - 2.5*alog10(1E-6*objectnew[num].irac2) ;
;   if objectnew[num].irac3 lt 0 then mag3 = -1. else mag3 = 5.16 - 2.5*alog10(1E-6*objectnew[num].irac3);
;   if objectnew[num].irac4 lt 0 then mag4 = -1. else mag4 = 4.52 - 2.5*alog10(1E-6*objectnew[num].irac4)
;   if objectnew[num].irac1 lt 0 then err1 = -1. else err1 =  0.05*(6.12 - 2.5*alog10(1E-6*objectnew[num].irac1))
;   if objectnew[num].irac2 lt 0 then err2 = -1. else err2 =  0.05*(5.64 - 2.5*alog10(1E-6*objectnew[num].irac2))
;   if objectnew[num].irac3 lt 0 then err3 = -1. else err3 =  0.05*(5.16 - 2.5*alog10(1E-6*objectnew[num].irac3))
;   if objectnew[num].irac4 lt 0 then err4 = -1. else err4 =  0.05*(4.52 - 2.5*alog10(1E-6*objectnew[num].irac4))

;reject some stars based on fwhm's
;reject the stars later
;   if objectnew[num].gfwhm gt 0 and objectnew[num].gfwhm lT 5.6 then begin
;      print, objectnew[num].gfwhm
      ;don't want these stars
;      starcount = starcount + 1
;   endif else begin
;      if       objectnew[num].umaga lt 0 and objectnew[num].umaga gt 28 and $
;        objectnew[num].gmaga lt 0 and objectnew[num].gmaga gt 28 and  $
;        objectnew[num].rmaga lt 0 and objectnew[num].rmaga gt 27.4 and $
;        objectnew[num].imaga lt 0 and objectnew[num].imaga gt 26.8 and $
;        objectnew[num].irac1 lt 0 and 8.296 - 2.5*alog10(1E-6*objectnew[num].irac1)  gt 21.7  and $
;        objectnew[num].irac2 lt 0 and 8.296 - 2.5*alog10(1E-6*objectnew[num].irac2)  gt 21.2 and $
;        objectnew[num].irac3 lt 0 and 8.296 - 2.5*alog10(1E-6*objectnew[num].irac3)  gt 20.2  and $
;        objectnew[num].irac4 lt 0 and 8.296 - 2.5*alog10(1E-6*objectnew[num].irac4)  gt 19.4  then begin

;   if objectnew[num].tmassj gt 0  and objectnew[num].tmassj lt 30 then tmassjab = -2.5*alog10(objectnew[num].tmassj) + 8.926 else tmassjab = objectnew[num].tmassj
;    if objectnew[num].tmassj gt 0  and objectnew[num].tmassh lt 30 then tmasshab = -2.5*alog10(objectnew[num].tmassh) + 8.926 else tmasshab = objectnew[num].tmassh
;   if objectnew[num].tmassj gt 0  and objectnew[num].tmassk lt 30 then tmasskab = -2.5*alog10(objectnew[num].tmassk) + 8.926 else tmasskab = objectnew[num].tmassk


 ;  if objectnew[num].irac1 lt 0 then mag1 = -1. else mag1 = 8.926 - 2.5*alog10(1E-6*objectnew[num].irac1) 
 ;  if objectnew[num].irac2 lt 0 then mag2 = -1. else mag2 = 8.926 - 2.5*alog10(1E-6*objectnew[num].irac2) 
 ;  if objectnew[num].irac3 lt 0 then mag3 = -1. else mag3 = 8.926- 2.5*alog10(1E-6*objectnew[num].irac3)
 ;  if objectnew[num].irac4 lt 0 then mag4 = -1. else mag4 = 8.926 - 2.5*alog10(1E-6*objectnew[num].irac4)
;(8.926 - 2.5*alog10(1E-6*objectnew[num].irac1)) -( 8.926 - 2.5*alog10(1E-6*objectnew[num].irac1 + 0.05*1E-6*objectnew[num].irac1))

;   if objectnew[num].irac1 lt 0 then err1 = -1. else err1 =  0.05*(8.926 - 2.5*alog10(1E-6*objectnew[num].irac1))
;   if objectnew[num].irac2 lt 0 then err2 = -1. else err2 =  0.05*(8.926 - 2.5*alog10(1E-6*objectnew[num].irac2));
;   if objectnew[num].irac3 lt 0 then err3 = -1. else err3 =  0.05*(8.926 - 2.5*alog10(1E-6*objectnew[num].irac3))
;   if objectnew[num].irac4 lt 0 then err4 = -1. else err4 =  0.05*(8.926 - 2.5*alog10(1E-6*objectnew[num].irac4))
;a = 1.4
;f = 1E-5
;err = 0.05*1E-5*a
;mag = 8.96 - 2.5*alog10(1E-5*a)
;magerr = (8.96 - 2.5*alog10(1E-5*a)) - (8.96 - 2.5*alog10(1E-5*a+0.05*1E-5*a))
;print, f, err, mag, magerr
;plot, objectnew.flux, objectnew.fluxerr, psym = 3, xrange=[0,1000], yrange=[0,50]

;newarr = objectnew.flux
;newarr2 = objectnew.flux
;j = 0
;for i = 0, n_elements(objectnew.flux) -1 do begin
;   if objectnew[i].fluxerr lt 0.05*objectnew[i].flux then begin
;      newarr[num]= objectnew[i].flux
;      newarr2[num] = objectnew[i].mips24
;      j = j + 1
;   endif;
;endfor
;newarr = newarr[0:j-1]
;newarr2 = newarr2[0:j-1]

;plot, newarr, newarr2, psym = 3, xrange = [0,1000], yrange=[0,1000]
;   if err1 gt 1000. then err1 = 1000.
;   if err2 gt 1000. then err2 = 1000.
;   if err3 gt 1000. then err3 = 1000.
;   if err4 gt 1000. then err4 = 1000.
