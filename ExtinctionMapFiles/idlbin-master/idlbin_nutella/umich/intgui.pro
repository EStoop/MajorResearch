PRO intgui
close, /all

;by hand pick a bunch of stars, here are there fwhm from
;SExtractor.cat
; find the average fwhm and sigma(fwhm)

fwhmarr = fltarr(15)
;A141
;fwhmarr=[3.65,3.71,3.51,3.44,3.39 ,3.32,3.34,3.34,3.45,3.64, 3.40,3.26,3.23,3.31,3.56,3.21,3.12,3.25,3.30,3.28,3.39]

;A3984
;fwhmarr=[4.14,4.06,3.91,4.06,3.96,3.85,3.96,3.77,3.86,4.04,3.85,4.42,4.05,4.07 ,4.08,3.93]

;A2556
fwhmarr=[3.29,3.29,3.24,3.18,3.16,3.18,3.30,3.48,3.30,3.40 ,3.46,3.58,3.56,3.49,3.43]
avfwhm = mean(fwhmarr)
sig = stddev(fwhmarr)
print, avfwhm, avfwhm+2*sig

;fits_read, "/n/Godiva2/jkrick/A141/A141.wcs.fits", data, header
;openr, lun,"/n/Godiva2/jkrick/A141/SExtractor.cat", /get_lun
;openw, gallun, "/n/Godiva2/jkrick/A141/galaxies1.cat",/get_lun
;openw, starlun, "/n/Godiva2/jkrick/A141/stars.cat",/get_lun

;fits_read, "/n/Godiva2/jkrick/A3984/A3984r.big.div.fits", data, header
;openr, lun,"/n/Godiva2/jkrick/A3984/SExtractor.cat", /get_lun
;openw, gallun, "/n/Godiva2/jkrick/A3984/galaxies3.cat",/get_lun
;openw, starlun, "/n/Godiva2/jkrick/A3984/stars.cat",/get_lun

fits_read, "/n/Godiva3/jkrick/A2556/A2556.large.fits", data, header
openr, lun,"/n/Godiva3/jkrick/A2556/SExtractor.cat", /get_lun
openw, gallun, "/n/Godiva3/jkrick/A2556/galaxies.cat",/get_lun
openw, starlun, "/n/Godiva3/jkrick/A2556/stars.cat",/get_lun

i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    readf, lun, num, x,y,a,b,e, flux, mag, magerr, isoarea, fwhm,theta, back

    isorad =  sqrt(isoarea/(!PI*(1-e)))    ;in pixels
    isorad = isorad *0.259                         ;in arcseconds
 
    musthave = 0.
   ;be careful of saturated stars
    ;A141;
;    IF (x GT 2920.) AND (x LT 3011.) AND ( y GT 2988.) AND ( y LT 3090.) THEN flux = -100.
;    IF (x GT 2116) AND (x LT 2160) AND ( y GT 2575) AND ( y LT 2603) THEN flux = -100
;    IF (x GT 2665) AND (x LT 2705) AND ( y GT 2648) AND ( y LT 2677) THEN flux = -100
;    IF (x GT 965) AND (x LT 999) AND ( y GT 2784) AND ( y LT 2813) THEN flux = -100
;    IF (x GT 1572) AND (x LT 1606) AND ( y GT 1067) AND ( y LT 1096) THEN flux = -100
;    IF (x GT 761) AND (x LT 790) AND ( y GT 767) AND ( y LT 790) THEN flux = -100
;    IF (x GT 642) AND (x LT 682) AND ( y GT 733) AND ( y LT 801) THEN flux = -100
;    IF (x GT 2091) AND (x LT 2179) AND ( y GT 572) AND ( y LT 668) THEN flux = -100
;    IF (x GT 980) AND (x LT 1013) AND ( y GT 3595) AND ( y LT 3627) THEN flux = -100
;   IF (x GT 564) AND (x LT 612) AND ( y GT 3107) AND ( y LT 3147) THEN flux = -100
;   IF (x GT 1924) AND (x LT 1948) AND ( y GT 2043) AND ( y LT 2051) THEN flux = -100
;   IF (x GT 1640) AND (x LT 1659) AND ( y GT 1863) AND ( y LT 1883) THEN flux = -100
;   IF (x GT 2225) AND (x LT 2245) AND ( y GT 1487) AND ( y LT 1507) THEN flux = -100
;   IF (x GT 1609) AND (x LT 1629) AND ( y GT 1022) AND ( y LT 1042) THEN flux = -100

;A3984
;    IF (x GT 588) AND (x LT 636) AND ( y GT 2291) AND ( y LT 2331) THEN flux = -100
;    IF (x GT 1020) AND (x LT 1060) AND ( y GT 2187) AND ( y LT 2267) THEN flux = -100
;    IF (x GT 788) AND (x LT 1020) AND ( y GT 684) AND ( y LT 948) THEN flux = -100
;    IF (x GT 2651) AND (x LT 2819) AND ( y GT 1190) AND ( y LT 1364) THEN flux = -100
;    IF (x GT 3019) AND (x LT 3067) AND ( y GT 2635) AND ( y LT 2707) THEN flux = -100
;    IF (x GT 2219) AND (x LT 2275) AND ( y GT 2867) AND ( y LT 2923) THEN flux = -100
;    IF (x GT 708) AND (x LT 772) AND ( y GT 724) AND ( y LT 780) THEN flux = -100
;    IF (x GT 470) AND (x LT 540) AND ( y GT 3531) AND ( y LT 3600) THEN flux = -100
;    IF (x GT 2364) AND (x LT 2384) AND ( y GT 2557) AND ( y LT 2577) THEN flux = -100
;    IF (x GT 1513) AND (x LT 1533) AND ( y GT 1922) AND ( y LT 1942) THEN flux = -100
;    IF (x GT 2887) AND (x LT 2907) AND ( y GT 631) AND ( y LT 651) THEN flux = -100
;    IF (x GT 2352) AND (x LT 2372) AND ( y GT 1566) AND ( y LT 1586) THEN flux = -100
;    IF (x GT 3148) AND (x LT 3168) AND ( y GT 2395) AND ( y LT 2415) THEN flux = -100
;    IF (x GT 3215) AND (x LT 3224) AND ( y GT 2661) AND ( y LT 2671) THEN flux = -100

;A2556
    IF (x GT 936.07) AND (x LT 976.07) AND (y GT 4026.82) AND (y LT 4066.82) THEN flux = -100
    IF (x GT 1766.35) AND (x LT 1806.35) AND (y GT 3692.39) AND (y LT 3732.39) THEN flux = -100
    IF (x GT 2613.47) AND (x LT 2653.47) AND (y GT 3666.58) AND (y LT 3706.58) THEN flux = -100
    IF (x GT 1944.28) AND (x LT 1984.28) AND (y GT 3524.14) AND (y LT 3564.14) THEN flux = -100
    IF (x GT 1963.51) AND (x LT 2003.51) AND (y GT 3418.11) AND (y LT 3458.11) THEN flux = -100
    IF (x GT 960.3) AND (x LT 1000.3) AND (y GT 3238.68) AND (y LT 3278.68) THEN flux = -100
    IF (x GT 214.18) AND (x LT 254.18) AND (y GT 2695.71) AND (y LT 2735.71) THEN flux = -100
    IF (x GT 616.45) AND (x LT 656.45) AND (y GT 2459.24) AND (y LT 2499.24) THEN flux = -100
    IF (x GT 1844.18) AND (x LT 1884.18) AND (y GT 2339.92) AND (y LT 2379.92) THEN flux = -100
    IF (x GT 2342.43) AND (x LT 2382.43) AND (y GT 2469.68) AND (y LT 2509.68) THEN flux = -100
    IF (x GT 2599.31) AND (x LT 2639.31) AND (y GT 2274.11) AND (y LT 2314.11) THEN flux = -100
   IF (x GT 1855.07) AND (x LT 1895.07) AND (y GT 2129.36) AND (y LT 2169.36) THEN flux = -100
    IF (x GT 751.29) AND (x LT 791.29) AND (y GT 1722.93) AND (y LT 1762.93) THEN flux = -100
    IF (x GT 89.93) AND (x LT 129.93) AND (y GT 1375.47) AND (y LT 1415.47) THEN flux = -100
    IF (x GT 109.19) AND (x LT 437.19) AND (y GT 816.9) AND (y LT 856.9) THEN flux = -100
    IF (x GT 320.04) AND (x LT 360.04) AND (y GT 445.31) AND (y LT 485.31) THEN flux = -100
    IF (x GT 651.1) AND (x LT 691.1) AND (y GT 476.03) AND (y LT 516.03) THEN flux = -100
    IF (x GT 2284.45) AND (x LT 2324.45) AND (y GT 430.42) AND (y LT 470.42) THEN flux = -100
    IF (x GT 977.46) AND (x LT 1017.46) AND (y GT 5.74) AND (y LT 45.74) THEN flux = -100
    IF (x GT 1104) AND (x LT 1144) AND (y GT 3600) AND (y LT 3630) THEN flux = -100
if (x GT  707.74 ) and (x LT  737.74 ) and (y gt  4018.28 ) and (y lt 4058.28 ) then flux = -100
if (x GT  1836.54 ) and (x LT  1866.54 ) and (y gt  4061.97 ) and (y lt 4101.97 ) then flux = -100
if (x GT  2632.28 ) and (x LT  2662.28 ) and (y gt  3339.22 ) and (y lt 3379.22 ) then flux = -100
if (x GT  2753.53 ) and (x LT  2783.53 ) and (y gt  3443.03 ) and (y lt 3483.03 ) then flux = -100
if (x GT  1246.78 ) and (x LT  1276.78 ) and (y gt  3273.09 ) and (y lt 3313.09 ) then flux = -100
if (x GT  1359.05 ) and (x LT  1389.05 ) and (y gt  3345.59 ) and (y lt 3385.59 ) then flux = -100
if (x GT  1967.94 ) and (x LT  1997.94 ) and (y gt  3126.99 ) and (y lt 3166.99 ) then flux = -100
if (x GT  1047.45 ) and (x LT  1077.45 ) and (y gt  2969.88 ) and (y lt 3009.88 ) then flux = -100
if (x GT  1216.01 ) and (x LT  1246.01 ) and (y gt  2920.68 ) and (y lt 2960.68 ) then flux = -100
if (x GT  2074.31 ) and (x LT  2104.31 ) and (y gt  2909.58 ) and (y lt 2949.58 ) then flux = -100
if (x GT  47.98 ) and (x LT  77.98 ) and (y gt  2679.39 ) and (y lt 2719.39 ) then flux = -100
if (x GT  1003.42 ) and (x LT  1033.42 ) and (y gt  2322.49 ) and (y lt 2362.49 ) then flux = -100
if (x GT  1123.48 ) and (x LT  1153.48 ) and (y gt  2242.09 ) and (y lt 2282.09 ) then flux = -100
if (x GT  795.46 ) and (x LT  825.46 ) and (y gt  947.44 ) and (y lt 987.44 ) then flux = -100
if (x GT  2188.66 ) and (x LT  2218.66 ) and (y gt  1365.74 ) and (y lt 1405.74 ) then flux = -100
if (x GT  1510.03 ) and (x LT  1540.03 ) and (y gt  666.54 ) and (y lt 706.54 ) then flux = -100

    IF (isoarea GT 6) AND (flux GT 0) AND (mag LT 25.5)THEN BEGIN ;is a real object

        IF ( e GE 0.1) OR (fwhm GT avfwhm+2*sig) THEN BEGIN ;is a galaxy
            IF (mag LT 20.3) AND (mag GE 14.3) AND (isorad GT 1.0)THEN begin
                xyad, header, x, y, ra, dec

                radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc            
                IF mag LE 17.0 THEN musthave = -1.5
                
                printf,gallun,      strcompress("@gvl"+string(i), /remove_all), "   ", $
                       strcompress(string(ihr)+":"+string(imin)+":"+string(xsec),/remove_all),   "     ",$
                       strcompress(string(ideg)+":"+string(imn)+":"+string(xsc),/remove_all), "   ",$
                                 musthave, fix(0), 1.0, 2, isorad + 4,isorad+4,0.0
;                printf, gallun, x, y
                i = i + 1
            ENDIF
        ENDIF

        IF ( e LT 0.1) AND (fwhm LT avfwhm + 2*sig) AND (mag LT 19) AND (mag GT 17)THEN BEGIN ;is a star
            xyad, header, x,y,ra,dec
            radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc
           printf,starlun,  strcompress("*s"+string(i), /remove_all), "   ", $
             strcompress(string(ihr)+":"+string(imin)+":"+string(xsec),/remove_all),   "     ",$
             strcompress(string(ideg)+":"+string(imn)+":"+string(xsc),/remove_all)
;            printf, starlun, x,y
            i = i + 1
        ENDIF

    ENDIF

    
ENDWHILE

close, lun
free_lun, lun
close,gallun
free_lun, gallun
close, starlun
free_lun, starlun
END
