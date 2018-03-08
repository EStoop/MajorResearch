pro clean_catalog

close, /all
restore, '/Users/jkrick/idlbin/object.sav'
openw, outlun, '/Users/jkrick/nep/singlelist', /get_lun

for j = 0l, n_elements(object.ra) - 1 do begin

;pick out all objects in the catalog with only one flux measurement.
   nflux = 0
   if object[j].umag gt 0 and object[j].umag ne 99 then nflux = nflux + 1
   if object[j].gmag gt 0 and object[j].gmag ne 99 then nflux = nflux + 1
   if object[j].rmag gt 0 and object[j].rmag ne 99 then nflux = nflux + 1
   if object[j].imag gt 0 and object[j].imag ne 99 then nflux = nflux + 1
   if object[j].acsmag gt 0 and object[j].acsmag ne 99 then nflux=nflux+1
   if object[j].flamjmag gt 0 and object[j].flamjmag ne 99 then nflux = nflux + 1
   if object[j].wircjmag gt 0 and object[j].wircjmag ne 99 then nflux = nflux + 1
   if object[j].wirchmag gt 0 and object[j].wirchmag ne 99 then nflux = nflux + 1
   if object[j].wirckmag gt 0 and object[j].wirckmag ne 99 then nflux = nflux + 1
   if object[j].tmassjmag gt 0 and object[j].tmassjmag ne 99 then nflux = nflux + 1
   if object[j].tmasshmag gt 0 and object[j].tmasshmag ne 99 then nflux = nflux + 1
   if object[j].tmasskmag gt 0 and object[j].tmasskmag ne 99 then nflux = nflux + 1
   if object[j].irac1mag gt 0 and object[j].irac1mag ne 99 then nflux = nflux + 1
   if object[j].irac2mag gt 0 and object[j].irac2mag ne 99 then nflux = nflux + 1
   if object[j].irac3mag gt 0 and object[j].irac3mag ne 99 then nflux = nflux + 1
   if object[j].irac4mag gt 0 and object[j].irac4mag ne 99 then nflux = nflux + 1
   if object[j].mips24mag gt 0 and object[j].mips24mag ne 99 then nflux = nflux + 1


   ;junk object
   if nflux lt 2 then begin
;      object[j].acsclassstar = 1
      printf, outlun, object[j].ra, object[j].dec
      printf, outlun, object[j].umag ,object[j].gmag , object[j].rmag , object[j].imag , object[j].acsmag, object[j].flamjmag ,object[j].wircjmag ,object[j].wirchmag ,object[j].wirckmag ,object[j].tmassjmag ,object[j].tmasshmag ,object[j].tmasskmag ,object[j].irac1mag ,object[j].irac2mag ,object[j].irac3mag ,object[j].irac4mag ,object[j].mips24mag 

   endif

;
endfor

;mark the stars
;have been careful to really only include stars and no galaxies
;I have therefore missed some stars and other unresolved objects, especially on the faint end,
;or regions where there is not as much multiwavelength data

targ0 = where(object.acsclassstar gt 0.98 or object.uclassstar gt 0.98 or object.gclassstar gt 0.98 or object.rclassstar gt 0.98 or object.iclassstar gt 0.98)
targ1 = where(object.acsfwhm lt 3.2 and object.acsellip lt 0.19 and object.acsfwhm gt 0 and object.acsellip gt 0)
;targ1 = where(object.iclassstar gt 0.90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19)
targ2 = where(object.iclassstar gt 0.90  and object.ifwhm lt 6.5 and object.iellip lt 0.20)
targ3 = where(object.ufluxmax gt 20 or object.gfluxmax gt 78 or object.rfluxmax gt 200. or object.ifluxmax gt 180 or object.acsfluxmax gt 124)
print, n_elements(targ0), n_elements(targ1)
targ = [targ0,targ1];,targ2,targ3]


px=[12.5,25, 25, 10, 12.5] ; acsmag 
py=[6e5 ,60, 35, 3e5, 6e5] ; ascisoarea 
py=alog10(py) ; switch to log space 

roi = Obj_New('IDLanROI', px, py) 

; test to see who is inside the polygon 
star=roi->ContainsPoints(object.acsmag,alog10(object.acsisoarea)) 

; at this point, inside the region is a star, i.e. star>1 
temp=where(star GT 0) 
targ = [targ0,temp];targ1];,targ2,targ3]

print, 'temp', n_elements(temp)


object[temp].acsclassstar = 1.0
;object[temp].uclassstar = 1.0
;object[temp].gclassstar = 1.0
;object[temp].rclassstar = 1.0
;object[temp].iclassstar = 1.0


close, outlun
free_lun, outlun

;save, object, filename= '/Users/jkrick/idlbin/object.sav'

end


;   if object[j].ufluxmax gt 20 then object[j].acsclassstar = 1
;   if object[j].gfluxmax gt 78 then object[j].acsclassstar = 1
;   if object[j].rfluxmax gt 200 then object[j].acsclassstar = 1
;   if object[j].ifluxmax gt 180 then object[j].acsclassstar = 1
  
;;   if object[j].acsfwhm lt 6.0  and object[j].acsfwhm gt 0. then object[j].acsclassstar = 1
;   if object[j].rfwhm lt 8.0 and  object[j].rfwhm gt 0. then object[j].acsclassstar = 1
