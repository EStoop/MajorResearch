PRO catalog_merge

close,/all

;instead could next time rad in the catalogs to an array, and then compare the arrays
;also to make a cmd, really don't want the stars in there.


i = 0
colorarr = fltarr(32000)
magarr = fltarr(32000)
 
openw, outlun, "/Users/jkrick/palomar/LFC/catalog/optical.cat", /get_lun
openr, lunr, "/Users/jkrick/palomar/LFC/catalog/SExtractor.r.cat", /get_lun
WHILE (NOT EOF(lunr)) DO BEGIN
    READF, lunr, o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magerrr, fwhmr, isoarear
    true = 0

    if magautor gt 0 and magautor lt 50 then begin

       openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat", /get_lun
       WHILE (NOT EOF(lung)) DO BEGIN
          readf, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag
          if (true EQ 0) then begin

        ;if the galaxies are the same, output the SExtractor parameters 
             IF (xwcsr LT xwcsg + 0.003) AND (xwcsr GT xwcsg - 0.0003) AND $
                (ywcsr LT ywcsg + 0.0003) AND (ywcsr GT ywcsg - 0.0003) and $
                (magautog gt 0) and (magautog LT 50) THEN BEGIN
                
                true = 1
                colorarr(i) = magautog - magautor
                magarr(i) =  magautor
                i = i + 1
                
             ENDIF 
          ENDIF

       ENDWHILE

       close,lung
       free_lun, lung

    ENDIF
       
  ENDWHILE
;meanarr = meanarr[0:i-1]

;print, mean(meanarr)

close,lunr



colorarr = colorarr(0:i-1)
magarr = magarr(0:i-1)

printf, outlun, colorarr
printf, outlun, magarr

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/palomar/LFC/catalog/g_r.ps", /portrait, xsize = 4, ysize = 4,/color

plot, magarr, colorarr, xthick=3, charthick=3, ythick=3,psym = 2

ps_close, /noprint, /noid
close, outlun
free_lun,outlun
free_lun,lunr

END



;V
;xshift = 623.88
;yshift = -656.7

;r
;xshift = 50.
;yshift = -20.

;	xshift	=	-210.9261 - (-124.4917)
;	yshift	=	-18.40577 - (953.8102)

;openr, lun1, "/n/Godiva6/jkrick/A3880/original/ccd3057.obsfile", /get_lun
;openw,lun1out, "/n/Godiva6/jkrick/A3880/original/ccd3057.obsfile2", /get_lun
;WHILE (NOT EOF(lun1)) DO BEGIN
;    readf, lun1, o,xcen,ycen,a, b, e,
;    flux_best,mag_best,isoarea,fwhm,pa,back,$
;    flux_aper,flux_isocor,flux_iso
;    readf, lun1, o, xcen, ycen, a , b, e, flux_iso, mag_iso, magerr, $
;      isoarea, fwhm, pa, back
;    true = 0;
;    print, "working on", xcen, ycen,fwhm
;    openr, lun2, "/n/Godiva6/jkrick/A3880/original/ccd9067.obsfile", /get_lun
;    WHILE (NOT EOF(lun2))DO BEGIN
;        readf, lun2, o2, xcen2, ycen2, a2 , b2, e2, flux_iso2, mag_iso2, $
;          magerr2, isoarea2, fwhm2, pa2, back2
;        IF (xcen2 - xshift GT xcen - 3.0) AND (xcen2 - xshift LT xcen+ 3.0) AND $  
;          (ycen2- yshift GT ycen- 3.0) AND (ycen2 - yshift LT ycen + 3.0) THEN BEGIN 
;;          AND   (e5 LT 0.4) AND (mag_iso5 LT 23) THEN BEGIN  ;376  1003
;           ;object is the same
;             true = 1
;             print, xcen2, ycen2, xcen, ycen
 ;            printf, lun1out,  fix(o2), xcen2, ycen2,xcen, ycen, mag_iso, magerr, 1.525641
;;                   magerr, isoarea, fwhm, pa, back
 ;         ENDIF 
 ;    ENDWHILE
 ;    close, lun2
 ;    free_lun, lun2
 ; ENDWHILE
;close,lun1
;close,lun1out
;free_lun, lun1
;free_lun, lun1out;


;END

;             dist = sqrt((1436. - xcenter)^2 + (1524.-ycenter)^2)
;             IF (dist GT 224) AND dist LT 1072 AND (24.6 -
;             2.5*alog10(isocorflux) LE 24.5) AND fwhm $
;             GT 5.1 AND isoarea GT 6 THEN BEGIN;


;             IF (24.3 - 2.5*alog10(isoco) GT 17.03) AND (24.3 - 2.5*alog10(isoco) LT 17.25) THEN begin


           ;object is a member of the cluster
;            meanarr[i] = (24.3 - 2.5*alog10(isoco)) - (24.6 - 2.5*alog10(isocorflux))
;                 printf,outlun,xcenter,ycenter, (24.3 - 2.5*alog10(isoco)),(24.6 - 2.5*alog10(isocorflux))
;                 meanarr[i] = 24.6 - 2.5*alog10(isocorflux)
;                 printf, outlun, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
            ;printf, outlun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
;            i = i + 1
;             ENDIF
