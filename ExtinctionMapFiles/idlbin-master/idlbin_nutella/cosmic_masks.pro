pro cosmic_masks

close, /all
;openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/bmasklist.txt', /get_lun  ;b = bright = < 21.5
;openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/dmasklist.txt', /get_lun ; d = dim = > 21.5
openw, outlun, '/Users/jkrick/conferences/aas2007/gallist.txt', /get_lun 

openw, outlun2, '/Users/jkrick/palomar/cosmic/slitmasks/slitobjmatch.txt',/get_lun

restore, '/Users/jkrick/idlbin/object.sav'
zphotarr = fltarr(10000)
eroarr = fltarr(10000)

hd=headfits('/Users/jkrick/hst/raw/wholeacs.fits')
rhd = headfits('/Users/jkrick/palomar/lfc/coadd_r.fits')
ihd = headfits('/Users/jkrick/palomar/lfc/coadd_i.fits')

priority = 0



   ;find the stars
   ;don't care if they are saturated

for i = 0l, n_elements(object.rmaga) -1 do begin

   if object[i].rmaga lt 18.5  and object[i].rmaga gt 0 and object[i].rfwhm lt 8.0  then begin
      radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
      printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
      printf, outlun2, priority, i
      priority = priority + 1
   endif

endfor

;---------------------------------
;find interesting galaxies
priority = 300


for i = 0l, n_elements(object.rmaga) -1 do begin

nflux = 0
if object[i].umaga gt 0 and object[i].umaga ne 99 then nflux = nflux + 1
if object[i].gmaga gt 0 and object[i].gmaga ne 99 then nflux = nflux + 1
if object[i].rmaga gt 0 and object[i].rmaga ne 99 then nflux = nflux + 1
if object[i].imaga gt 0 and object[i].imaga ne 99 then nflux = nflux + 1
if object[i].acsmag gt 0 and object[i].acsmag ne 99 then nflux=nflux+1
if object[i].flamjmag gt 0 and object[i].flamjmag ne 99 then nflux = nflux + 1
if object[i].wircjmag gt 0 and object[i].wircjmag ne 99 then nflux = nflux + 1
if object[i].wirchmag gt 0 and object[i].wirchmag ne 99 then nflux = nflux + 1
if object[i].wirckmag gt 0 and object[i].wirckmag ne 99 then nflux = nflux + 1
if object[i].tmassjmag gt 0 and object[i].tmassjmag ne 99 then nflux = nflux + 1
if object[i].tmasshmag gt 0 and object[i].tmasshmag ne 99 then nflux = nflux + 1
if object[i].tmasskmag gt 0 and object[i].tmasskmag ne 99 then nflux = nflux + 1
if object[i].irac1mag gt 0 and object[i].irac1mag ne 99 then nflux = nflux + 1
if object[i].irac2mag gt 0 and object[i].irac2mag ne 99 then nflux = nflux + 1
if object[i].irac3mag gt 0 and object[i].irac3mag ne 99 then nflux = nflux + 1
if object[i].irac4mag gt 0 and object[i].irac4mag ne 99 then nflux = nflux + 1
if object[i].mips24mag gt 0 and object[i].mips24mag ne 99 then nflux = nflux + 1


   if object[i].rmaga le 28.0 and object[i].acsclassstar lt 1. then begin ; and object[i].zphot gt 0.05 and object[i].zphot le 1.2 then begin

      ; red things
      if object[i].irac1mag gt 0. and object[i].rmaga - object[i].irac1mag ge 4.0 then begin    ;4.0 is real definition
         adxy,rhd, object[i].ra, object[i].dec, x, y
         print, "   ero", i, object[i].rmaga, object[i].irac1mag, object[i].zphot
;         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
;         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
;         printf, outlun2, priority, i
         priority = priority + 1
      endif

      ; distant things
;      if object[i].wircjmag lt 99 and object[i].wircjmag - object[i].wirckmag ge 1.6 then begin   ;2.3 is real definition
;           adxy,hd, object[i].ra, object[i].dec, x, y
;       print,x, y,  "    DRG", object[i].rmaga, object[i].wircjmag
;      endif

      if object[i].zphot gt 1.0 and object[i].prob gt 80  and nflux gt 2 and object[i].imaga lt 21.5 then begin
         adxy,ihd, object[i].ra, object[i].dec, x, y
        print,"distant" , i, object[i].ifwhm, object[i].imaga
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority, i

         priority = priority + 1

;         print, x, y, "   high z" , object[i].zphot, object[i].prob       ;want to check these sed's
      endif

      ;70 micron
      if object[i].mips70mag gt 0 then begin
         adxy,hd, object[i].ra, object[i].dec, x, y
         print,"   s70micron", i, object[i].mips70mag
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority, i
         priority = priority + 1
      endif


      ;blue things
      if object[i].irac1mag lt 99 and object[i].irac1mag gt 0 and object[i].rmaga gt 0  and object[i].irac1mag - object[i].rmaga gt 1.5 then begin
         adxy,hd, object[i].ra, object[i].dec, x, y
         print,  "   blue", i,  object[i].rmaga
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority, i
         priority = priority + 1
         
      endif

   endif
endfor

;---------------------------------

   ; find the random galaxies
priority = 1000
priority2 = 2000
priority3 = 3000
priority4 = 4000
priority5 = 5000
count = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0

j = 0
zarr = fltarr(30000)
for i = 0l, n_elements(object.rmaga) -1 do begin

nflux = 0
if object[i].umaga gt 0 and object[i].umaga ne 99 then nflux = nflux + 1
if object[i].gmaga gt 0 and object[i].gmaga ne 99 then nflux = nflux + 1
if object[i].rmaga gt 0 and object[i].rmaga ne 99 then nflux = nflux + 1
if object[i].imaga gt 0 and object[i].imaga ne 99 then nflux = nflux + 1
if object[i].acsmag gt 0 and object[i].acsmag ne 99 then nflux=nflux+1
if object[i].flamjmag gt 0 and object[i].flamjmag ne 99 then nflux = nflux + 1
if object[i].wircjmag gt 0 and object[i].wircjmag ne 99 then nflux = nflux + 1
if object[i].wirchmag gt 0 and object[i].wirchmag ne 99 then nflux = nflux + 1
if object[i].wirckmag gt 0 and object[i].wirckmag ne 99 then nflux = nflux + 1
if object[i].tmassjmag gt 0 and object[i].tmassjmag ne 99 then nflux = nflux + 1
if object[i].tmasshmag gt 0 and object[i].tmasshmag ne 99 then nflux = nflux + 1
if object[i].tmasskmag gt 0 and object[i].tmasskmag ne 99 then nflux = nflux + 1
if object[i].irac1mag gt 0 and object[i].irac1mag ne 99 then nflux = nflux + 1
if object[i].irac2mag gt 0 and object[i].irac2mag ne 99 then nflux = nflux + 1
if object[i].irac3mag gt 0 and object[i].irac3mag ne 99 then nflux = nflux + 1
if object[i].irac4mag gt 0 and object[i].irac4mag ne 99 then nflux = nflux + 1
if object[i].mips24mag gt 0 and object[i].mips24mag ne 99 then nflux = nflux + 1


   if object[i].rmaga le 22.0 and object[i].rmaga gt 0 and object[i].acsclassstar lt 1. then begin ; and object[i].zphot gt 0.05 and object[i].zphot le 1.2 then begin
      zarr[j] = object[i].zphot
      j = j + 1
      if object[i].zphot gt 0.05 and object[i].zphot lt 0.2 then begin
         if object[i].prob gt 0 or nflux ge 5 then begin
            radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
            priority = priority + 1
            printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
            printf, outlun2, priority, i
               
            count = count  + 1
         endif
      endif
        if object[i].zphot gt 0.2 and object[i].zphot lt 0.4 then begin
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         priority2 = priority2 + 1
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',priority2, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority2, i
               
         count2 = count2  + 1
      endif
      
        if object[i].zphot gt 0.4 and object[i].zphot lt 0.6 then begin
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         priority3 = priority3 + 1
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',priority3, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority3, i
               
         count3 = count3  + 1
      endif
      
        if object[i].zphot gt 0.6 and object[i].zphot lt 0.8 then begin
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         priority4 = priority4 + 1
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',priority4, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority4, i
      
         count4 = count4  + 1
      endif
      
        if object[i].zphot gt 0.8 and object[i].zphot lt 1.2 then begin
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
         priority5 = priority5 + 1
         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',priority5, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
         printf, outlun2, priority5, i

         count5 = count5  + 1
      endif
      

   endif

endfor
print, "there are ", count, " objects with mr < 22 and z < 0.2"
print, "there are ", count2, " objects with mr < 22 and 0.2 < z > 0.4"
print, "there are ", count3, " objects with mr < 22 and 0.4 < z > 0.6"
print, "there are ", count4, " objects with mr < 22 and 0.6 < z > 0.8"
print, "there are ", count5, " objects with mr < 22 and 0.8 < z > 1.2"

zarr = zarr[0:j-1]
plothist,  zarr, xhist, yhist, bin = 0.05, /noprint, xrange=[0,2]
;--------------------------------------------------
;add the targeted galaxies to the list
priority = 200
;read in Jason's list.
listjason =[    679 ,       980] 
for i = 0, n_elements(listjason) - 1 do begin
   radec, object[listjason(i)].ra, object[listjason(i)].dec, ihr, imin, xsec, ideg, imn, xsc
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[listjason(i)].rmaga
   printf, outlun2, priority, listjason(i)

   priority = priority + 1
endfor
;ashby's list
;70 micron by hand with r mag < 22
list70 = [      21274 ,      1366  ,       1480 ,       1391 ,     10847  ,      1234  ,      1235  ,        949 ,        475 ,       563  ,       535    ,       561  ]
for i = 0, n_elements(list70) - 1 do begin
   print, "i 70", list70(i)
   radec, object[list70(i)].ra, object[list70(i)].dec, ihr, imin, xsec, ideg, imn, xsc
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[list70(i)].rmaga
   printf, outlun2, priority, list70(i)

   priority = priority + 1
endfor

;possible cluster galaxy. only one bright enough
listc = [      10790 ]
   radec, object[listc].ra, object[listc].dec, ihr, imin, xsec, ideg, imn, xsc
   printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[listc].rmaga
   printf, outlun2, priority, listc

   priority = priority + 1


close, outlun
free_lun, outlun
close, outlun2
free_lun, outlun2

end

;make a list of 250 random galaxies 

;how many objects are detected in r band and have mr < 22 and have a redshift in the range 0.05 to 1.2? 
; ignoring acs data for now
;need to be optically detected because this is an optical spectrograph
;do I want to add a requirement that the phot z be accurate?  I think it is ok to have some slop.
;will need to weed out the stars


;what is the redshift distribution of those objects?
;plothist, zphotarr, xhist, yhist, bin=0.1, /noprint, xrange=[0.1,2]

; how do I choose which ones in the bins where there are way more than 20?
; randomly?
;those which have detections in the most other bands?
;do I need to worry about the selection function?

;make a list of ra, dec, rmag. 
;also make a list of x,y to overlay on rband image for optical confirmation
;do I need stars to orient?

;-----------------------------------------
;need a list of 50 interesting targets

;ero's = r - 3.6 >=4.0; none




