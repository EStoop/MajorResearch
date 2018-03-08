pro addz, zobject, object

print, 'starting addz'
;close, /all

;!P.multi = [0,2,2]

;numoobjects = 74000
;colorarr = fltarr(numoobjects)
;magarr = fltarr(numoobjects)

;zobject = replicate({zob, zra:0D,zdec:0D,zfluxauto:0D,zmagauto:0D,zmagerrauto:0D,zfluxaper:0D,zmagaper:0D,zmagerraper:0D,zfluxbest:0D,zmagbest:0D,zmagerrbest:0D,zfwhm:0D,zisoarea:0D, zellip:0D,zfluxmax:0D, zclassstar:0D, zflags:0D},numoobjects)

;restore, '/Users/jkrick/idlbin/object.sav'
;----------------------------------------------------------------

;p = 0.
;openr, luni, "/Users/jkrick/mmt/SExtractor.cat", /get_lun

;fits_read, '/Users/jkrick/mmt/IRACCF.z.coadd.cov.fits', zdata, zhead

;WHILE (NOT EOF(luni)) DO BEGIN
;   READF, luni, NUMBER,X_WORLD      ,Y_WORLD      ,X_IMAGE,Y_IMAGE,FLUX_AUTO,MAG_AUTO,MAGERR_AUTO,FLUX_APER,MAG_APER,MAGERR_APER,FLUX_BEST,MAG_BEST,MAGERR_BEST,FWHM_IMAGE,ISOAREA_IMAGE,FLUX_MAX,ELLIPTICITY,CLASS_STAR,FLAGS
;   if zdata[x_image, y_image] gt 1 then begin
;      zobject[p] ={zob, x_world, y_world, flux_auto, MAG_AUTO,MAGERR_AUTO,FLUX_APER,MAG_APER,MAGERR_APER,FLUX_BEST,MAG_BEST,MAGERR_BEST,FWHM_IMAGE,ISOAREA_IMAGE,ELLIPTICITY,flux_max,CLASS_STAR,FLAGS}
;      p = p + 1
 ;  endif

; endwhile
;zobject =zobject[0:p - 1]
;print, "there are ",p,"z objects"
;close, luni
;free_lun, luni
;----------------------------------------------------------------



;----------------------------------------------------------
;match z to object


; create initial arrays
m=n_elements(zobject.zra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999

print,'Matching z to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
   dist=sphdist( zobject[q].zra, zobject[q].zdec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
;   print, zobject[q].zra, zobject[q].zdec, sep, object[ind].ra, object[ind].dec

   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
   endif 
endfor

print,"Finished at "+systime()
nonmatched = (where(mmatch lt 0))
print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))
openw, outlunz, '/Users/jkrick/mmt/znonmatch.reg', /get_lun
printf, outlunz, 'fk5'
for num2=0, n_elements(nonmatched) - 1 do begin
   
printf, outlunz, 'circle (', zobject(nonmatched[num2]).zra, zobject(nonmatched[num2]).zdec,' 10" )'


endfor

close, outlunz
free_lun, outlunz


;undefine, zobject
end


