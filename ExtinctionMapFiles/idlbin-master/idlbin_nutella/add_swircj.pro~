pro addswircj

;have a shallow wide catalog and a deep small area catalog
;first match the wide catalog to populate all of the objects with
;their J mag.
;second match the deep catalog to replace shallow with deep where
;possible.



restore, '/Users/jkrick/idlbin/objectnew.sav'

readcol, '/Users/jkrick/mmt/jband/NEPswmos.txt', ra, dec, xcenter, ycenter,magaper, magapererr, magauto, magautoerr, kronrad, flags, fwhm, elong, ellip, classstar, format='A'



swircjobject = {swircjxcenter:0D, swircjycenter:0D, swircjra:0D, swircjdec:0D, swircjmagaper:0D, swircjmagapererr:0D, swircjmagauto:0D, swircjmagautoerr:0D,swircjkronrad:0D, swircjflags:0D, swircjfwhm:0D, swircjelong:0D, swircjellip:0D , swircjclassstar:0D, swircjnndist:0D, swircjmatchdist:0D, swircjmatch:0D}

b = replicate(swircjobject, n_elements(objectnew.ra))
objectnew=struct_addtags(objectnew, b)

objectnew.swircjmatch = -999
objectnew.swircjmatchdist = 999





;do some matching
m = n_elements(ra)
mmatch = fltarr(m)
mmatch[*] = -999

for q=0,m -1 do begin
   dist=sphdist( ra(q), dec(q),objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin 
      if objectnew[ind].swircjmatch ge 0 then begin   ;if there was a previous match
         if sep lt objectnew[ind].swircjmatchdist then begin   ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].swircjmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].swircjra = ra(q)
            objectnew[ind].swircjdec = dec(q)
            objectnew[ind].swircjxcenter=xcenter(q)
            objectnew[ind].swircjycenter=ycenter(q)
            objectnew[ind].swircjmagaper=magaper(q)
            objectnew[ind].swircjmagapererr=magapererr(q)
            objectnew[ind].swircjmagauto=magauto(q)
            objectnew[ind].swircjmagautoerr=magautoerr(q)
            objectnew[ind].swircjkronrad=kronrad(q)
            objectnew[ind].swircjflags=flags(q)
            objectnew[ind].swircjfwhm=fwhm(q)
            objectnew[ind].swircjelong=elong(q)
            objectnew[ind].swircjellip=ellip(q)
            objectnew[ind].swircjclassstar=classstar(q)
            objectnew[ind].swircjnndist=-99
            objectnew[ind].xmatch = q
            objectnew[ind].xmatchdist = sep

         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
            objectnew[ind].swircjra = ra(q)
            objectnew[ind].swircjdec = dec(q)
            objectnew[ind].swircjxcenter=xcenter(q)
            objectnew[ind].swircjycenter=ycenter(q)
            objectnew[ind].swircjmagaper=magaper(q)
            objectnew[ind].swircjmagapererr=magapererr(q)
            objectnew[ind].swircjmagauto=magauto(q)
            objectnew[ind].swircjmagautoerr=magautoerr(q)
            objectnew[ind].swircjkronrad=kronrad(q)
            objectnew[ind].swircjflags=flags(q)
            objectnew[ind].swircjfwhm=fwhm(q)
            objectnew[ind].swircjelong=elong(q)
            objectnew[ind].swircjellip=ellip(q)
            objectnew[ind].swircjclassstar=classstar(q)
            objectnew[ind].swircjnndist=-99
            objectnew[ind].xmatch = q
            objectnew[ind].xmatchdist = sep

      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

;print out any that do not have a match

openw, outlunred, '/Users/jkrick/nep/highz/robertson.reg', /get_lun
printf, outlunred, 'fk5'
for num= 0,  n_elements(nonmatched) -1. do printf,outlunred,  'circle( ', ra(nonmatched[num]), dec(nonmatched[num]),  ' 3")', magaper(nonmatched[num])
close, outlunred
free_lun, outlunred



save, objectnew, filename='/Users/jkrick/idlbin/objectnew_swirc.sav'


end
