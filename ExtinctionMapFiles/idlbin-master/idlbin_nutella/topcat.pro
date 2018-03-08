pro topcat

restore, '/Users/jkrick/idlbin/object.sav'
openw, outlun,  '/Users/jkrick/nep/catalog.txt', /get_lun

count = 1l
starcount = 0

printf, outlun, '#num  ra  dec  umaga   gmaga  rmaga  imaga   acsmag  zmagbest  jmagab  wircjmagab  tmassjmag   wirchmagab  tmasshmag  wirckmagab  tmasskmag  irac1mag irac2mag irac3mag irac4mag    umagerra  gmagerra  rmagerra  imagerra  acsmagerr  zmagerrbest  flamjmagerr  wircjmagerr tmassjerr  wirchmagerr  tmassherr  wirckmagerr tmasskerr err1 err2 err3 err4'

for num = 0l, n_elements(object.uxcenter) -1 do begin
;AB magnitudes
   ;don't convert things which are already -1, or -99
   
   if object[num].flamjmag gt 0 and object[num].flamjmag ne 99 then begin
      fab = 1594.*10^(object[num].flamjmag/(-2.5))
      jmagab = -2.5*alog10(fab) +8.926
   endif else begin
      jmagab = object[num].flamjmag
   endelse

   if object[num].wircjmag gt 0 and object[num].wircjmag ne 99 then begin
      wircjab = 1594.*10^(object[num].wircjmag/(-2.5))
      wircjmagab = -2.5*alog10(wircjab) +8.926
   endif else begin
      wircjmagab = object[num].wircjmag
   endelse

   if object[num].wirchmag gt 0 and object[num].wirchmag ne 99 then begin
      wirchab = 1024.*10^(object[num].wirchmag/(-2.5))
      wirchmagab = -2.5*alog10(wirchab) +8.926
   endif else begin
      wirchmagab = object[num].wirchmag
   endelse

   if object[num].wirckmag gt 0 and object[num].wirckmag ne 99 then begin
      wirckab = 666.8*10^(object[num].wirckmag/(-2.5))
      wirckmagab = -2.5*alog10(wirckab) +8.926
   endif else begin
      wirckmagab = object[num].wirckmag
   endelse

   if object[num].irac1flux lt 0 then err1 = -1. else err1 = 0.05
   if object[num].irac2flux lt 0 then err2 = -1. else err2 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[num].irac2)
   if object[num].irac3flux lt 0 then err3 = -1. else err3 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[num].irac3)
   if object[num].irac4flux lt 0 then err4 = -1. else err4 =  0.05;8.926 - 2.5*alog10(0.05E-6*object[num].irac4)
   
   if object[num].imagerra gt 1000. then object[num].imagerra = 1000.
   if object[num].gmagerra gt 1000. then object[num].gmagerra = 1000.
   if object[num].rmagerra gt 1000. then object[num].rmagerra = 1000.
   if object[num].umagerra gt 1000. then object[num].umagerra = 1000.
   if object[num].zmagerrbest gt 1000. then object[num].zmagerrbest = 1000.
   
   if object[num].tmassjmag lt 0 then tmassjerr = 99 else tmassjerr = 0.02
   if object[num].tmasshmag lt 0 then tmassherr = 99 else tmassherr = 0.02
   if object[num].tmasskmag lt 0 then tmasskerr = 99 else tmasskerr = 0.02

   nflux = 0
   if object[num].umaga gt 0 and object[num].umaga ne 99 then nflux = nflux + 1
   if object[num].gmaga gt 0 and object[num].gmaga ne 99 then nflux = nflux + 1
   if object[num].rmaga gt 0 and object[num].rmaga ne 99 then nflux = nflux + 1
   if object[num].imaga gt 0 and object[num].imaga ne 99 then nflux = nflux + 1
   if object[num].acsmag gt 0 and object[num].acsmag ne 99 then nflux=nflux+1
   if object[num].zmagbest gt 0 and object[num].zmagbest ne 99 then nflux = nflux +1
   if object[num].flamjmag gt 0 and object[num].flamjmag ne 99 then nflux = nflux + 1
   if object[num].wircjmag gt 0 and object[num].wircjmag ne 99 then nflux = nflux + 1
   if object[num].wirchmag gt 0 and object[num].wirchmag ne 99 then nflux = nflux + 1
   if object[num].wirckmag gt 0 and object[num].wirckmag ne 99 then nflux = nflux + 1
   if object[num].tmassjmag gt 0 and object[num].tmassjmag ne 99 then nflux = nflux + 1
   if object[num].tmasshmag gt 0 and object[num].tmasshmag ne 99 then nflux = nflux + 1
   if object[num].tmasskmag gt 0 and object[num].tmasskmag ne 99 then nflux = nflux + 1
   if object[num].irac1mag gt 0 and object[num].irac1mag ne 99 then nflux = nflux + 1
   if object[num].irac2mag gt 0 and object[num].irac2mag ne 99 then nflux = nflux + 1
   if object[num].irac3mag gt 0 and object[num].irac3mag ne 99 then nflux = nflux + 1
   if object[num].irac4mag gt 0 and object[num].irac4mag ne 99 then nflux = nflux + 1
   if object[num].mips24mag gt 0 and object[num].mips24mag ne 99 then nflux = nflux + 1

   if nflux gt 0 then  printf, outlun, format='(I10,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
                 num, object[num].ra, object[num].dec, object[num].umaga , object[num].gmaga, object[num].rmaga, $
                 object[num].imaga,  object[num].acsmag, object[num].zmagbest, jmagab, wircjmagab, $
                 object[num].tmassjmag,  wirchmagab, object[num].tmasshmag, wirckmagab, $
                 object[num].tmasskmag, object[num].irac1mag,object[num].irac2mag,$
                 object[num].irac3mag,object[num].irac4mag,   $
                 object[num].umagerra, object[num].gmagerra, $
                 object[num].rmagerra, object[num].imagerra, object[num].acsmagerr, object[num].zmagerrbest, object[num].flamjmagerr, $
                 object[num].wircjmagerr,tmassjerr, object[num].wirchmagerr, tmassherr, object[num].wirckmagerr,$
                 tmasskerr,err1,err2,err3,err4
;took out acs until photometry gets sorted out.

         count = count + 1
  

endfor
print, "count", count
print, "num", num
close, outlun
free_lun, outlun

end

