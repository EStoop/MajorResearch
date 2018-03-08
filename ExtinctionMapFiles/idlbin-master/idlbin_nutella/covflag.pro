;-------------------------------------------------------------------------------------------
; takes a list of input images and a list of coords, checks for non-NaN
;  coverage in the images and outputs a flag
PRO COVFLAG
close,/all
imagelist = '/Users/jkrick/nep/covlist'
;positionlist = '/Users/jkrick/nep/testposition'
restore, '/Users/jkrick/idlbin/object.sav'

readcol,imagelist,images, format="a"
nfiles=n_elements(images)
print, nfiles
;readcol,positionlist,cntr,ra,dec,format="a"
ra = object.ra
dec = object.dec
cntr = findgen(n_elements(object.ra))
npos=n_elements(cntr)

coverage=intarr(npos)
coverage[*]=0

; start looping through the input images
for i=0l,nfiles-1 do begin

   print,images[i]
; for speed, just read the header
   if i eq 4 then bcd_hd=headfits(images[i],exten=2) else bcd_hd=headfits(images[i])
   naxis1=sxpar(bcd_hd,'NAXIS1')-2
   naxis2=sxpar(bcd_hd,'NAXIS2')-2
   crval1=sxpar(bcd_hd,'CRVAL1')
   crval2=sxpar(bcd_hd,'CRVAL2')
   print, naxis1, naxis2
; get the x and y positions of all points relative to current image
   adxy,bcd_hd,ra,dec,x,y
   x=round(x)
   y=round(y)
   
; check to see if they are on image
   onimage=where(((x  LT naxis1) AND (x  GT 0) AND (y  LT naxis2) AND (y  GT 0)),count)
;   onimage=where(((x +20. LT naxis1) AND (x - 20. GT 2) AND (y + 20. LT naxis2) AND (y - 20. GT 2)),count)
;   print,"count", count

; begin the loop for checking coverage
   if count NE 0 then begin
      
; read the actual image
      if i eq 4 then bcd=readfits(images[i],bcd_hd,exten_no=2) else bcd=readfits(images[i],bcd_hd)

      for j=0l,count-1 do begin
         state=0
         coverage(onimage(j)) = 0

         if  finite (bcd(x(onimage(j)),y(onimage(j))) ) GT 0. and (bcd(x(onimage(j)),y(onimage(j))) ne 0.0)  then   state=state+1

;         if  finite (bcd(x(onimage(j)),y(onimage(j))) ) GT 0. and (bcd(x(onimage(j)),y(onimage(j))) ne 0.0)  then   state=state+1
;         if  finite (bcd(20+x(onimage(j)),y(onimage(j)))) gt 0. and (bcd(20+x(onimage(j)),y(onimage(j))) ne 0.0) then state=state+1
;         if finite (bcd(-20+x(onimage(j)),y(onimage(j))) ) gt 0. and  (bcd(-20+x(onimage(j)),y(onimage(j))) ne 0.) then state=state+1
;         if finite (bcd(x(onimage(j)),20+y(onimage(j))) ) gt 0. and (bcd(x(onimage(j)),20+y(onimage(j))) ne 0.) then state=state+1
;         if finite (bcd(x(onimage(j)),-20+y(onimage(j))) ) gt 0. and  (bcd(x(onimage(j)),-20+y(onimage(j))) ne 0.) then state=state+1

         if (state eq 1) then begin
;         if (state GT 4) then begin
            coverage(onimage(j))=1
         
         endif
;         if i eq 4 then print, "coverage", ra(onimage(j)), dec(onimage(j)), coverage(onimage(j))

; end the loop over overlapping objects
      endfor
; end the branch if there are overlaps on image
   endif 
; end the loop over images

   ;try to get the ones which are off the image to have coverage=0
   for lm =long(0), n_elements(object.ra) - 1 do begin
      if  x(lm) lt 0. or y(lm) lt 0. or x(lm) gt naxis1 - 20. or y(lm) gt naxis1 - 20. then  coverage(lm) =0
   endfor

   for k=0l,npos-1 do begin
      if i eq 0 and  coverage(k) lt 1. then object[(k)].umag = -1
      if i eq 1 and  coverage(k) lt 1. then object[(k)].gmag = -1
      if i eq 2 AND coverage(k) lt 1. then object[(k)].rmag = -1
      if i eq 3 AND coverage(k) lt 1. then object[(k)].imag = -1
      if i eq 4 and coverage(k) lt 1. then object[(k)].acsmag = -1
      if i eq 5 AND coverage(k) lt 1. then object[(k)].flamjmag = -1
      if i eq 6 AND coverage(k) lt 1. then object[(k)].wircjmag = -1
      if i eq 7 AND coverage(k) lt 1. then object[(k)].wirchmag = -1
      if i eq 8 AND coverage(k) lt 1. then object[(k)].wirckmag = -1
      if i eq 9 AND coverage(k) lt 1. then begin
         object[(k)].irac1mag = -1
         object[(k)].irac1flux = -1
      endif
      if i eq 10 AND coverage(k) lt 1. then begin
         object[(k)].irac2mag = -1
         object[(k)].irac2flux = -1
      endif
      if i eq 11 AND coverage(k) lt 1. then begin
         object[(k)].irac3mag = -1
         object[(k)].irac3flux = -1
      endif
      if i eq 12 AND coverage(k) lt 1. then begin
         object[(k)].irac4mag = -1
         object[(k)].irac4flux = -1
      endif
      if i eq 13 AND coverage(k) lt 1. then begin
         object[(k)].mips24mag = -1
         object[(k)].mips24flux = -1
      endif
      if i eq 14 and coverage(k) lt 1. then begin
         object[(k)].zmagbest = -1
         object[(k)].zmagauto = -1
         object[(k)].zmagaper = -1
      endif

   endfor

endfor

; dump the output
openw,unit,'/Users/jkrick/nep/test.txt',/get_lun

for i=0l,npos-1 do begin
   printf,unit,cntr(i),ra(i),dec(i),coverage(i)
;   if coverage(i) lt 1. then object[cntr(i)].gmaga = -1
endfor

free_lun,unit


;object(where(object.umag lt -90)).umag = 99.0
;object(where(object.gmag lt -90)).gmag = 99.0
;object(where(object.rmag lt -90)).rmag = 99.0
;object(where(object.imag lt -90)).imag = 99.0
;object(where(object.zmagbest lt -90)).zmagbest = 99.0
object(where(object.irac1mag lt -90)).irac1mag = 99.0
object(where(object.irac2mag lt -90)).irac2mag = 99.0
object(where(object.irac3mag lt -90)).irac3mag = 99.0
object(where(object.irac4mag lt -90)).irac4mag = 99.0
;object(where(object.mips24mag lt -90)).mips24mag = 99.0

save, object, filename= '/Users/jkrick/idlbin/object.test.sav'


END

;if the object is off the frame
;but all off frame object should have coverag=0 anyway, since coverage array starts out filled with zeros
;   for lm =0, n_elements(object.ra) - 1 do begin

;      if x lt 0. or y lt 0. or x gt naxis1 - 10. or y gt naxis1 - 10. then  coverage(lm) =0


;  
;      if  (bcd(x(onimage(j)),y(onimage(j))) ne 0.0)  then   state=state+1
;      if   (bcd(1+x(onimage(j)),y(onimage(j))) ne 0.0) then state=state+1
;      if  (bcd(-1+x(onimage(j)),y(onimage(j))) ne 0.) then state=state+1
;      if (bcd(x(onimage(j)),1+y(onimage(j))) ne 0.) then state=state+1
;      if  (bcd(x(onimage(j)),-1+y(onimage(j))) ne 0.) then state=state+1
;         print, "state", state
