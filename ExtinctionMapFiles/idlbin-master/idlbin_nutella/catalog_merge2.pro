PRO catalog_merge

close,/all


numoobjects = 32000
colorarr = fltarr(numoobjects)
magarr = fltarr(numoobjects)

rgalaxy = replicate({object, xcenter:0D,ycenter:0D,fluxa:0D,maga:0D,magerra:0D,fwhm:0D,isoarea:0D},numoobjects)
ggalaxy = replicate({gobject, gxcenter:0D,gycenter:0D,gfluxa:0D,gmaga:0D,gmagerrb:0D,gfwhm:0D,gisoarea:0D},numoobjects)
galaxy= replicate({what, rra:0D, rdec:0D, g_r:0D, rmag:0D, gmag:0D}, numoobjects)

openw, outlun, "/Users/jkrick/palomar/LFC/catalog/optical.cat", /get_lun
;fluxmax, fwhm  limits for saturated stars, star/galaxy separation
;r = 150, 5.6
;i = 110  ??
;g=76, 8.4
;u = 19??
i = 0

openr, lunr, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat", /get_lun
WHILE (NOT EOF(lunr)) DO BEGIN
    READF, lunr, o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magerrr, fwhmr, isoarear, fluxmaxr

    if (magautor LT 50) and (fwhmr gt 5.6) and  (fluxmaxr lt 76) then begin
       rgalaxy[i] ={object, xwcsr, ywcsr,fluxautor, magautor, magerrr,fwhmr, isoarear}
       i = i + 1
;       printf, outlun, xcenterr, ycenterr
    ENDIF

 endwhile
rgalaxy = rgalaxy[0:i - 1]
print, "there are i,",i," r galaxies"
close, lunr
free_lun, lunr


l = 0
m = 0.
openr, lung, "/Users/jkrick/palomar/LFC/catalog/SExtractor.u.cat", /get_lun
WHILE (NOT EOF(lung)) DO BEGIN
    READF, lung, o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag, fluxmaxg

    if (fluxautog gt 0 ) and (fwhmg gt 8.4) and (fluxmaxg lt 20) then begin
       ggalaxy[l] ={object, xwcsg, ywcsg,fluxautog, magautog, magerrg,fwhmg, isoareag}
;       print, "wprking on galaxy", l, xcenterg, ycenterg, xwcsg, ywcsg

       FOR j = 0, i-1, 1 do begin
          IF (rgalaxy[j].xcenter LT xwcsg + 0.0003) AND (rgalaxy[j].xcenter GT xwcsg - 0.0003) AND $
             (rgalaxy[j].ycenter  LT ywcsg + 0.0003) AND (rgalaxy[j].ycenter GT ywcsg - 0.0003) then begin
             
             galaxy[m] = {what, rgalaxy[j].xcenter, rgalaxy[j].ycenter, magautog - rgalaxy[j].maga, rgalaxy[j].maga, magautog}
             m = m + 1.
             printf,outlun, xcenterg, ycenterg, rgalaxy[j].xcenter, rgalaxy[j].ycenter, magautog, fwhmg
          ENDIF

       ENDFOR

       l = l + 1
    ENDIF

 endwhile
ggalaxy = ggalaxy[0:l - 1]
print, "there are i,",l," g galaxies"
close, lung
free_lun, lung

galaxy = galaxy[0:m-1]
print, "there are i,",m," matched galaxies"


!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/palomar/LFC/catalog/u_g.ps", /portrait, xsize = 4, ysize = 4,/color

plot, galaxy.rmag, galaxy.g_r, xthick=3, charthick=3, ythick=3,psym = 2, xrange=[10,30],yrange=[-5,5],$
ytitle="u-g", xtitle="g"

ps_close, /noprint, /noid

close, outlun
free_lun,outlun
;free_lun,lunr


undefine, rgalaxy
undefine, ggalaxy
undefine, galaxy
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
