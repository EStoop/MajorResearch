pro addwirck
restore, '/Users/jkrick/idlbin/objectnew.sav'
;;;;;;;

;;;;;;;;


numoobjects = 10000

wirckobject = replicate({wirckob, wirckxcenter:0D, wirckycenter:0D, wirckra:0D, wirckdec:0D, wirckmag:0D, wirckmagerr:0D, wirckfwhm:0D, wirckisoarea:0D, wirckellip:0D , wircknndist:0D, wirckmatchdist:0D, wirckmatch:0D},numoobjects)

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/wirc_2008/k/kband_cat_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wirckobject[p] ={wirckob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
wirckobject = wirckobject[0:p - 1]
print, "there are ",p," wirc h objects"
close, luni
free_lun, luni



;determine nndist
for nn = 0l, n_elements(wirckobject.wirckra) - 1 do begin
   a = where(sqrt( (wirckobject[nn].wirckra - wirckobject.wirckra)^2 + (wirckobject[nn].wirckdec - wirckobject.wirckdec)^2 ) le 0.0008 )
   wirckobject[nn].wircknndist = n_elements(a)
endfor


;add new columns to catalog to accomodate old photometry.

wirckobjectold = replicate({wirckobold, wirckxcenter_old:0D, wirckycenter_old:0D, wirckra_old:0D, wirckdec_old:0D, wirckmag_old:0D, wirckmagerr_old:0D, wirckfwhm_old:0D, wirckisoarea_old:0D, wirckellip_old:0D , wircknndist_old:0D, wirckmatchdist_old:0D, wirckmatch_old:0D},n_elements(objectnew.ra))

objectnew = struct_addtags(objectnew, wirckobjectold)

objectnew.wirckra_old = objectnew.wirckra
objectnew.wirckdec_old = objectnew.wirckdec
objectnew.wirckxcenter_old = objectnew.wirckxcenter
objectnew.wirckycenter_old = objectnew.wirckycenter
objectnew.wirckmag_old = objectnew.wirckmag
objectnew.wirckmagerr_old = objectnew.wirckmagerr
objectnew.wirckfwhm_old = objectnew.wirckfwhm
objectnew.wirckisoarea_old = objectnew.wirckisoarea
objectnew.wirckellip_old = objectnew.wirckellip
objectnew.wircknndist_old = objectnew.wircknndist
objectnew.wirckmatch_old = objectnew.wirckmatch
objectnew.wirckmatchdist_old = objectnew.wirckmatchdist


;then clear the matches 
objectnew.wirckra = 0
objectnew.wirckdec = 0
objectnew.wirckxcenter = 0
objectnew.wirckycenter = 0
objectnew.wirckmag = 99.
objectnew.wirckmagerr = 99.
objectnew.wirckfwhm = 0
objectnew.wirckisoarea = 0
objectnew.wirckellip = 0
objectnew.wircknndist = 0
objectnew.wirckmatch = -999
objectnew.wirckmatchdist = 999


;ready to do some matching
;considered converting from vega to AB but plythyperz etc. already
;expect vega.
m=n_elements(wirckobject.wirckra)
mmatch=fltarr(m)
mmatch[*]=-999
print,'Matching wirck to object'
print,"Starting at "+systime()
dist=mmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( wirckobject[q].wirckra, wirckobject[q].wirckdec,objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if objectnew[ind].wirckmatch GE 0 then begin           ;if there was a previous match
         if sep lt objectnew[ind].wirckmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].wirckmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].wirckra = wirckobject[q].wirckra
            objectnew[ind].wirckdec = wirckobject[q].wirckdec
            objectnew[ind].wirckxcenter = wirckobject[q].wirckxcenter
            objectnew[ind].wirckycenter = wirckobject[q].wirckycenter
            objectnew[ind].wirckmag = wirckobject[q].wirckmag
            objectnew[ind].wirckmagerr = wirckobject[q].wirckmagerr
            objectnew[ind].wirckfwhm = wirckobject[q].wirckfwhm
            objectnew[ind].wirckisoarea = wirckobject[q].wirckisoarea
            objectnew[ind].wirckellip = wirckobject[q].wirckellip
            objectnew[ind].wircknndist =wirckobject[q].wircknndist
            objectnew[ind].wirckmatch = q
            objectnew[ind].wirckmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].wirckra = wirckobject[q].wirckra
         objectnew[ind].wirckdec = wirckobject[q].wirckdec
         objectnew[ind].wirckxcenter = wirckobject[q].wirckxcenter
         objectnew[ind].wirckycenter = wirckobject[q].wirckycenter
         objectnew[ind].wirckmag = wirckobject[q].wirckmag
         objectnew[ind].wirckmagerr = wirckobject[q].wirckmagerr
         objectnew[ind].wirckfwhm = wirckobject[q].wirckfwhm
         objectnew[ind].wirckisoarea = wirckobject[q].wirckisoarea
         objectnew[ind].wirckellip = wirckobject[q].wirckellip
         objectnew[ind].wircknndist =wirckobject[q].wircknndist
         objectnew[ind].wirckmatchdist = sep
         objectnew[ind].wirckmatch = q
      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
print,"Finished at "+systime()


;print out any that do not have a match
openw, outlunred, '/Users/jkrick/palomar/wirc/wirc_2008/h/nonmatch.reg', /get_lun
printf, outlunred, 'fk5'

for num= 0,  n_elements(nonmatched) -1. do begin
   
   printf, outlunred, 'circle( ', wirckobject[nonmatched[num]].wirckra,wirckobject[nonmatched[num]].wirckdec, ' 3")'
endfor
close, outlunred
free_lun, outlunred


;diagnostics

plot, objectnew.wirckmag_old, objectnew.wirckmag, psym = 2, thick = 3, xtitle='old', ytitle='new', xrange=[12,22], yrange=[12,22], xstyle = 1, ystyle = 1
oplot, findgen(30), findgen(30)



save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'


end
