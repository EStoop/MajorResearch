pro addwirch
restore, '/Users/jkrick/idlbin/objectnew.sav'
;;;;;;;
;consider fixing object 5915 which has Nan's while I am
;re-writing the catalog
;;;;;;;;


numoobjects = 10000

wirchobject = replicate({wirchob, wirchxcenter:0D, wirchycenter:0D, wirchra:0D, wirchdec:0D, wirchmag:0D, wirchmagerr:0D, wirchfwhm:0D, wirchisoarea:0D, wirchellip:0D , wirchnndist:0D, wirchmatchdist:0D, wirchmatch:0D},numoobjects)

p = 0
openr, luni, "/Users/jkrick/palomar/wirc/wirc_2008/h/hband_cat_nohead.txt", /get_lun
WHILE (NOT EOF(luni)) DO BEGIN
   READF, luni, junk,x_image, y_image,junk,junk,junk,junk,junk,junk,mag_best, magerr_best, junk,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk 

   wirchobject[p] ={wirchob, x_image, y_image, ra, dec, mag_best, magerr_best, fwhm_image, isoareaimage, ell,0.,0.,-1.}
   p = p + 1
   
 endwhile
wirchobject = wirchobject[0:p - 1]
print, "there are ",p," wirc h objects"
close, luni
free_lun, luni



;determine nndist
for nn = 0l, n_elements(wirchobject.wirchra) - 1 do begin
   a = where(sqrt( (wirchobject[nn].wirchra - wirchobject.wirchra)^2 + (wirchobject[nn].wirchdec - wirchobject.wirchdec)^2 ) le 0.0008 )
   wirchobject[nn].wirchnndist = n_elements(a)
endfor


;add new columns to catalog to accomodate old photometry.

wirchobjectold = replicate({wirchobold, wirchxcenter_old:0D, wirchycenter_old:0D, wirchra_old:0D, wirchdec_old:0D, wirchmag_old:0D, wirchmagerr_old:0D, wirchfwhm_old:0D, wirchisoarea_old:0D, wirchellip_old:0D , wirchnndist_old:0D, wirchmatchdist_old:0D, wirchmatch_old:0D},n_elements(objectnew.ra))

objectnew = struct_addtags(objectnew, wirchobjectold)

objectnew.wirchra_old = objectnew.wirchra
objectnew.wirchdec_old = objectnew.wirchdec
objectnew.wirchxcenter_old = objectnew.wirchxcenter
objectnew.wirchycenter_old = objectnew.wirchycenter
objectnew.wirchmag_old = objectnew.wirchmag
objectnew.wirchmagerr_old = objectnew.wirchmagerr
objectnew.wirchfwhm_old = objectnew.wirchfwhm
objectnew.wirchisoarea_old = objectnew.wirchisoarea
objectnew.wirchellip_old = objectnew.wirchellip
objectnew.wirchnndist_old = objectnew.wirchnndist
objectnew.wirchmatch_old = objectnew.wirchmatch
objectnew.wirchmatchdist_old = objectnew.wirchmatchdist


;then clear the matches 
objectnew.wirchra = 0
objectnew.wirchdec = 0
objectnew.wirchxcenter = 0
objectnew.wirchycenter = 0
objectnew.wirchmag = 99.
objectnew.wirchmagerr = 99.
objectnew.wirchfwhm = 0
objectnew.wirchisoarea = 0
objectnew.wirchellip = 0
objectnew.wirchnndist = 0
objectnew.wirchmatch = -999
objectnew.wirchmatchdist = 999


;ready to do some matching
;considered converting from vega to AB but plythyperz etc. already
;expect vega.
m=n_elements(wirchobject.wirchra)
mmatch=fltarr(m)
mmatch[*]=-999
print,'Matching wirch to object'
print,"Starting at "+systime()
dist=mmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( wirchobject[q].wirchra, wirchobject[q].wirchdec,objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      if objectnew[ind].wirchmatch GE 0 then begin           ;if there was a previous match
         if sep lt objectnew[ind].wirchmatchdist then begin  ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].wirchmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].wirchra = wirchobject[q].wirchra
            objectnew[ind].wirchdec = wirchobject[q].wirchdec
            objectnew[ind].wirchxcenter = wirchobject[q].wirchxcenter
            objectnew[ind].wirchycenter = wirchobject[q].wirchycenter
            objectnew[ind].wirchmag = wirchobject[q].wirchmag
            objectnew[ind].wirchmagerr = wirchobject[q].wirchmagerr
            objectnew[ind].wirchfwhm = wirchobject[q].wirchfwhm
            objectnew[ind].wirchisoarea = wirchobject[q].wirchisoarea
            objectnew[ind].wirchellip = wirchobject[q].wirchellip
            objectnew[ind].wirchnndist =wirchobject[q].wirchnndist
            objectnew[ind].wirchmatch = q
            objectnew[ind].wirchmatchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].wirchra = wirchobject[q].wirchra
         objectnew[ind].wirchdec = wirchobject[q].wirchdec
         objectnew[ind].wirchxcenter = wirchobject[q].wirchxcenter
         objectnew[ind].wirchycenter = wirchobject[q].wirchycenter
         objectnew[ind].wirchmag = wirchobject[q].wirchmag
         objectnew[ind].wirchmagerr = wirchobject[q].wirchmagerr
         objectnew[ind].wirchfwhm = wirchobject[q].wirchfwhm
         objectnew[ind].wirchisoarea = wirchobject[q].wirchisoarea
         objectnew[ind].wirchellip = wirchobject[q].wirchellip
         objectnew[ind].wirchnndist =wirchobject[q].wirchnndist
         objectnew[ind].wirchmatchdist = sep
         objectnew[ind].wirchmatch = q
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
   
   printf, outlunred, 'circle( ', wirchobject[nonmatched[num]].wirchra,wirchobject[nonmatched[num]].wirchdec, ' 3")'
endfor
close, outlunred
free_lun, outlunred


;diagnostics

plot, objectnew.wirchmag_old, objectnew.wirchmag, psym = 2, thick = 3, xtitle='old', ytitle='new', xrange=[12,22], yrange=[12,22], xstyle = 1, ystyle = 1
oplot, findgen(30), findgen(30)



save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'


end
