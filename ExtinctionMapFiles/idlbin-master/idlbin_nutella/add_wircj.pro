pro addwircj

;need to change this back to objectnew.sav

restore, '/Users/jkrick/idlbin/objectnew.sav'
;restore, '/Users/jkrick/idlbin/objectnew.test.sav'
;;;;;;;
;consider fixing object 5915 which has Nan's while I am
;re-writing the catalog
;if it is still a problem, they are in irac4fluxerr and irac4magerr
;;;;;;;;


numoobjects = 10000

wircjobject = replicate({wircjob, wircjxcenter:0D, wircjycenter:0D, wircjra:0D, wircjdec:0D, wircjmag:0D, wircjmagerr:0D, wircjfwhm:0D, wircjisoarea:0D, wircjellip:0D , wircjnndist:0D, wircjmatchdist:0D, wircjmatch:0D},numoobjects)

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/wirc_2008/j/jband_cat_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wircjobject[p] ={wircjob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
wircjobject = wircjobject[0:p - 1]
print, "there are ",p," wirc j objects"
close, luni
free_lun, luni



;determine nndist
for nn = 0l, n_elements(wircjobject.wircjra) - 1 do begin
   a = where(sqrt( (wircjobject[nn].wircjra - wircjobject.wircjra)^2 + (wircjobject[nn].wircjdec - wircjobject.wircjdec)^2 ) le 0.0008 )
   wircjobject[nn].wircjnndist = n_elements(a)
endfor


;add new columns to catalog to accomodate old photometry.

wircjobjectold = replicate({wircjobold, wircjxcenter_old:0D, wircjycenter_old:0D, wircjra_old:0D, wircjdec_old:0D, wircjmag_old:0D, wircjmagerr_old:0D, wircjfwhm_old:0D, wircjisoarea_old:0D, wircjellip_old:0D , wircjnndist_old:0D, wircjmatchdist_old:0D, wircjmatch_old:0D},n_elements(objectnew.ra))

objectnew = struct_addtags(objectnew, wircjobjectold)

objectnew.wircjra_old = objectnew.wircjra
objectnew.wircjdec_old = objectnew.wircjdec
objectnew.wircjxcenter_old = objectnew.wircjxcenter
objectnew.wircjycenter_old = objectnew.wircjycenter
objectnew.wircjmag_old = objectnew.wircjmag
objectnew.wircjmagerr_old = objectnew.wircjmagerr
objectnew.wircjfwhm_old = objectnew.wircjfwhm
objectnew.wircjisoarea_old = objectnew.wircjisoarea
objectnew.wircjellip_old = objectnew.wircjellip
objectnew.wircjnndist_old = objectnew.wircjnndist
objectnew.wircjmatch_old = objectnew.wircjmatch
objectnew.wircjmatchdist_old = objectnew.wircjmatchdist


;then clear the matches 
objectnew.wircjra = 0
objectnew.wircjdec = 0
objectnew.wircjxcenter = 0
objectnew.wircjycenter = 0
objectnew.wircjmag = 99.
objectnew.wircjmagerr = 99.
objectnew.wircjfwhm = 0
objectnew.wircjisoarea = 0
objectnew.wircjellip = 0
objectnew.wircjnndist = 0
objectnew.wircjmatch = -999
objectnew.wircjmatchdist = 999


;ready to do some matching
;considered converting from vega to AB but plythyperz etc. already
;expect vega.
m=n_elements(wircjobject.wircjra)
mmatch=fltarr(m)
mmatch[*]=-999
print,'Matching wircj to object'
print,"Starting at "+systime()
dist=mmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( wircjobject[q].wircjra, wircjobject[q].wircjdec,objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if objectnew[ind].wircjmatch GE 0 then begin           ;if there was a previous match
         if sep lt objectnew[ind].wircjmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].wircjmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].wircjra = wircjobject[q].wircjra
            objectnew[ind].wircjdec = wircjobject[q].wircjdec
            objectnew[ind].wircjxcenter = wircjobject[q].wircjxcenter
            objectnew[ind].wircjycenter = wircjobject[q].wircjycenter
            objectnew[ind].wircjmag = wircjobject[q].wircjmag
            objectnew[ind].wircjmagerr = wircjobject[q].wircjmagerr
            objectnew[ind].wircjfwhm = wircjobject[q].wircjfwhm
            objectnew[ind].wircjisoarea = wircjobject[q].wircjisoarea
            objectnew[ind].wircjellip = wircjobject[q].wircjellip
            objectnew[ind].wircjnndist =wircjobject[q].wircjnndist
            objectnew[ind].wircjmatch = q
            objectnew[ind].wircjmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].wircjra = wircjobject[q].wircjra
         objectnew[ind].wircjdec = wircjobject[q].wircjdec
         objectnew[ind].wircjxcenter = wircjobject[q].wircjxcenter
         objectnew[ind].wircjycenter = wircjobject[q].wircjycenter
         objectnew[ind].wircjmag = wircjobject[q].wircjmag
         objectnew[ind].wircjmagerr = wircjobject[q].wircjmagerr
         objectnew[ind].wircjfwhm = wircjobject[q].wircjfwhm
         objectnew[ind].wircjisoarea = wircjobject[q].wircjisoarea
         objectnew[ind].wircjellip = wircjobject[q].wircjellip
         objectnew[ind].wircjnndist =wircjobject[q].wircjnndist
         objectnew[ind].wircjmatchdist = sep
         objectnew[ind].wircjmatch = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()


;print out any that do not have a match
openw, outlunred, '/Users/jkrick/palomar/wirc/wirc_2008/j/nonmatch.reg', /get_lun
printf, outlunred, 'fk5'

for num= 0,  n_elements(nonmatched) -1. do begin
   
   printf, outlunred, 'circle( ', wircjobject[nonmatched[num]].wircjra,wircjobject[nonmatched[num]].wircjdec, ' 3")'
endfor
close, outlunred
free_lun, outlunred


;diagnostics
ps_open, filename='/Users/jkrick/palomar/wirc/wirc_2008/j/oldvsnew.ps',/portrait,/square,/color
plot, objectnew.wircjmag_old, objectnew.wircjmag, psym = 2, thick = 3, xtitle='old', ytitle='new', xrange=[12,22], yrange=[12,22], xstyle = 1, ystyle = 1
oplot, findgen(30), findgen(30)


ps_close, /noprint,/noid
save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'



end
