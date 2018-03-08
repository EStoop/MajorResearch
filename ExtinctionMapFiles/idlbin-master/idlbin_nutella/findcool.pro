pro findcool

close, /all
;openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/bmasklist.txt', /get_lun  ;b = bright = < 21.5
;openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/dmasklist.txt', /get_lun ; d = dim = > 21.5
restore, '/Users/jkrick/idlbin/object.sav'

hd=headfits('/Users/jkrick/hst/raw/wholeacs.fits')
rhd = headfits('/Users/jkrick/palomar/lfc/coadd_r.fits')
ihd = headfits('/Users/jkrick/palomar/lfc/coadd_i.fits')

 
;---------------------------------
;find interesting galaxies
for i = 0l, 3001 -1 do begin

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


   if object[i].rmaga le 28.0 and object[i].acsclassstar lt 1. and nflux gt 5 and object[i].acsmag lt 25 and object[i].acsmag gt 0 then begin ; and object[i].zphot gt 0.05 and object[i].zphot le 1.2 then begin

      ; red things
      if object[i].irac1mag gt 0. and object[i].rmaga - object[i].irac1mag ge 5.0 then begin    ;4.0 is real definition
         adxy,rhd, object[i].ra, object[i].dec, x, y
         print, "   ero", i, object[i].rmaga, object[i].irac1mag, object[i].zphot
;         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
;         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
 ;        priority = priority + 1
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
;         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
;         priority = priority + 1

;         print, x, y, "   high z" , object[i].zphot, object[i].prob       ;want to check these sed's
      endif

      ;70 micron
      if object[i].mips70mag gt 0 then begin
         adxy,hd, object[i].ra, object[i].dec, x, y
;         print,"   s70micron", i, object[i].mips70mag
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
;         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
;;         priority = priority + 1
      endif


      ;blue things
      if object[i].irac1mag lt 99 and object[i].irac1mag gt 0 and object[i].rmaga gt 0  and object[i].irac1mag - object[i].rmaga gt 1.5 then begin
         adxy,hd, object[i].ra, object[i].dec, x, y
         print,  "   blue", i,  object[i].rmaga
         radec, object[i].ra, object[i].dec, ihr, imin, xsec, ideg, imn, xsc
;         printf, outlun,  format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', priority, ihr, imin, xsec, ideg, imn, xsc, object[i].rmaga
;         priority = priority + 1
         
      endif


      if object[i].spt eq 1 and object[i].prob gt 80 then print, "arp220", i, object[i].ra, object[i].dec
     if object[i].spt eq 6 and object[i].prob gt 80 then print, "mrk231", i, object[i].ra, object[i].dec
     if object[i].spt eq 7 and object[i].prob gt 80 then print, "QS) 1", i, object[i].ra, object[i].dec
     if object[i].spt eq 8 and object[i].prob gt 80 then print, "QSO 2", i, object[i].ra, object[i].dec
     if object[i].spt eq 14 and object[i].prob gt 80 then print, "sey 2", i, object[i].ra, object[i].dec
     if object[i].spt eq 15 and object[i].prob gt 80 then print, "torus", i, object[i].ra, object[i].dec
     if object[i].spt eq 9 and object[i].prob gt 80 then print, "s0", i, object[i].ra, object[i].dec


   endif
endfor

;---------------------------------

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




