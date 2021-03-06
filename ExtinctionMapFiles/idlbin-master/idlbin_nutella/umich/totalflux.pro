PRO totalflux

close,/all

openr, lun,"/n/Godiva1/jkrick/A4059/original/elliptab", /get_lun

sum= 0
sumc= 0
area0 = 0.D
area = 0.D
areac0 = 0.D
areac = 0.D
intens=0.D
i = 0
rvir = (25.9 *223.)
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - areac0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        fluxc0 = areac*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF

    IF area GT 0 THEN begin
        sum = sum + area*intens
        sumc = sumc + areac*intens
        area0 = !Pi * sma^2*(1.-ellip)
        areac0 = !Pi * sma^2    ;circle test
    ENDIF
    i = i + 1
    ellipf = ellip   ;final ellipticity
    smaf = sma       ;final sma
    areaf = area0    ;final ellipses area
;    print, sma, intens, ellip, area, sum
ENDWHILE
print, smaf,areaf, ellipf
sum = sum - flux0
sumc = sumc - fluxc0

close, lun
free_lun, lun


OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2rtab', /GET_LUN

;read in the radial profile
rows= 20
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr1 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr1(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

start = [0.03,200.0]

;exponential fitting
result = MPFITFUN('exponential',radius[0:j-1],counts[0:j-1], Verr1[0:j-1], start)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100.
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

;add to the sum from points outside of the end radius already summed
area0 = areaf
bins = 100.
FOR n = 1, bins , 1 DO BEGIN
    sma = (rvir- smaf)/bins * n + smaf
    area = !Pi * sma^2*(1.-ellipf) - areaf
    rad = sma - (rvir-smaf)/(2*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    areaf = area + areaf
;    print, "A4059",sma, area, intens, sum
ENDFOR


m = ((22.04 - 2.5*alog10(sum)) - 36.64) 
l = 100^((4.63-m)/5)
;mc = ((22.04 - 2.5*alog10(sumc)) - 36.64) 
;lc = 100^((4.63-mc)/5)
;-----------------
;r galaxies
openr, lun,"/n/Godiva1/jkrick/A4059/original/ellipgalrtab", /get_lun

sum= 0
area0 = 0.D
areac0 = 0.D
areac = 0.D
area = 0.D
intens=0.D
i = 0
rvir = (34.9 *126.)

WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
    ellipf = ellip   ;final ellipticity
    smaf = sma       ;final sma
    areaf = area0    ;final ellipses area

ENDWHILE
close, lun
free_lun, lun

sum = sum - flux0

OPENR, lun,'/n/Godiva1/jkrick/A4059/original/gal2rtab', /GET_LUN

;read in the radial profile
rows= 20
rradius = FLTARR(rows)
rcounts = FLTARR(rows)
rerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius(j) = r
      rcounts(j) = c 
      rerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun

result = MPFITFUN('exponential',rradius,rcounts, err, start, /quiet)   

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


;add to the sum from points outside of the end radius already summed
area0 = areaf
bins = 100.
FOR n = 1, bins , 1 DO BEGIN
    sma = (rvir- smaf)/bins * n + smaf
    area = !Pi * sma^2*(1.-ellipf) - areaf
    rad = sma - (rvir-smaf)/(2*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    areaf = area + areaf
    print, "A4059",sma, area, intens, sum
ENDFOR


mg = ((24.6 - 2.5*alog10(sum)) - 39.71) 
lg = 100^((4.63-mg)/5)


;---------------------------------------------------------------------

openr, lun,"/n/Godiva1/jkrick/A4059/original/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0


close, lun
free_lun, lun



OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2Btab', /GET_LUN

;read in the radial profile
rows= 20 - 2
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]


;exponential fitting
result = MPFITFUN('exponential',radiusB,countsB, VerrB, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

m = ((22.19 - 2.5*alog10(sum)) - 36.64) 
lB = 100^((5.47-m)/5)

;---------
print, "cluster, VICL, Vgalaxies, rICL, rgalaxies"

print, "A4059", lB , 100^((5.47-(((22.19 - 2.5*alog10(5367)) - 36.64)))/5),lg,$
l, 100^((4.63-(((22.04 - 2.5*alog10(23949.2)) - 36.64)))/5)
;print, " A4059", l, lB, lc, l/lc

;-----------------------------------------------------------------
openr, lun,"/n/Godiva6/jkrick/A3880/original/elliptab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
areac0 = 0.D
areac = 0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN begin
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF

    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0
close, lun
free_lun, lun

OPENR, lun,'/n/Godiva6/jkrick/A3880/original/icl2rtab', /GET_LUN

;read in the radial profile
rows= 17 -2
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun


;exponential fitting
result = MPFITFUN('exponential',radius[0:j-1],counts[0:j-1], Verr, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.04 - 2.5*alog10(sum)) - 37.07) 
l = 100^((4.63-m)/5)

;-----------
openr, lun,"/n/Godiva6/jkrick/A3880/original/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF

    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0
close, lun
free_lun, lun

OPENR, lun,'/n/Godiva6/jkrick/A3880/original/icl2Btab', /GET_LUN

;read in the radial profile
rows= 10 
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]

;exponential fitting
result = MPFITFUN('exponential',radiusB,countsB, VerrB, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.19 - 2.5*alog10(sum)) - 37.07) 
lB = 100^((5.47-m)/5)

;-------------
print, "A3880", lB , 100^((5.47-(((22.19 - 2.5*alog10(3344.5)) - 37.07)))/5),$
l, 100^((4.63-(((22.04 - 2.5*alog10(11645.7)) - 37.07)))/5)

;print, "A3880", l,lB

;-----------------------------------------------------------------

openr, lun,"/n/Godiva4/jkrick/A2734/original/elliptab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
areac0 = 0.D
areac = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF

    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0
close, lun
free_lun, lun

OPENR, lun,'/n/Godiva4/jkrick/A2734/original/icl2rtab', /GET_LUN

;read in the radial profile
rows= 17 -2
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

;exponential fitting
result = MPFITFUN('exponential',radius[0:j-1],counts[0:j-1], Verr, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.04 - 2.5*alog10(sum)) - 37.22) 
l = 100^((4.63-m)/5)

;-----------------
openr, lun,"/n/Godiva4/jkrick/A2734/original/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF

    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva4/jkrick/A2734/original/icl2Btab', /GET_LUN

;read in the radial profile
rows= 12 -2
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]

;exponential fitting
result = MPFITFUN('exponential',radiusB,countsB, VerrB, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.19 - 2.5*alog10(sum)) - 37.22) 
lB = 100^((5.47-m)/5)
;--------
print, "A2734", lB , 100^((5.47-(((22.19 - 2.5*alog10(2631.2)) - 37.22)))/5),$
l, 100^((4.63-(((22.04 - 2.5*alog10(14063.9)) - 37.22)))/5)
;print, "A2734", l, lB

;-----------------------------------------------------------------

openr, lun,"/n/Godiva3/jkrick/A2556/original/elliptab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
areac0 = 0.D
areac = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva3/jkrick/A2556/original/icl2rtab', /GET_LUN

;read in the radial profile
rows= 15 - 3
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr1 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr1(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

;exponential fitting
result = MPFITFUN('exponential',radius[0:j-3],counts[0:j-3], Verr[0:j-3], start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.04 - 2.5*alog10(sum)) - 38.0) 
l = 100^((4.63-m)/5)


;--------------

openr, lun,"/n/Godiva3/jkrick/A2556/original/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva3/jkrick/A2556/original/icl2Btab', /GET_LUN

;read in the radial profile
rows= 10
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]
VerrB = VerrB[0:j-1]

;exponential fitting
start1=[0.0002, 36000]
result = MPFITFUN('exponential',radiusB,countsB, VerrB, start1,/quiet)    ;ICL
;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.19 - 2.5*alog10(sum)) - 38.0) 
lB = 100^((5.47-m)/5)
;----------
print, "A2556", lB , 100^((5.47-(((22.19 - 2.5*alog10(1421.8)) - 38.0)))/5),$
l, 100^((4.63-(((22.04 - 2.5*alog10(7232.)) - 38.0)))/5)
;print, "A2556", l,lB
;-----------------------------------------------------------------

openr, lun,"/n/Godiva4/jkrick/A4010/original/elliptab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
areac0 = 0.D
areac = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun


OPENR, lun,'/n/Godiva4/jkrick/A4010/original/icl2rtab', /GET_LUN

;read in the radial profile
rows= 29 -3
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun


;exponential fitting
result = MPFITFUN('exponential',radius[0:j-1],counts[0:j-1], Verr[0:J-1], start, /quiet)    ;ICL


;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((22.04 - 2.5*alog10(sum)) - 38.21) 
l = 100^((4.63-m)/5)

;---------------------------


openr, lun,"/n/Godiva4/jkrick/A4010/original/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva4/jkrick/A4010/original/icl2Btab', /GET_LUN

;read in the radial profile
rows= 32 -8
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]

;exponential fitting
result = MPFITFUN('exponential',radiusB,countsB, VerrB, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

m = ((22.19 - 2.5*alog10(sum)) - 38.21) 
lB = 100^((5.47-m)/5)

;---
print, "A4010", lB , 100^((5.47-(((22.19 - 2.5*alog10(1242.5)) - 38.21)))/5),$
l, 100^((4.63-(((22.04 - 2.5*alog10(5824.1)) - 38.21)))/5)
;print, "A4010", l,lB
;-----------------------------------------------------------------

openr, lun,"/n/Godiva1/jkrick/A3888/final/elliptab", /get_lun

sum= 0
area0 = 0.D
areac0 = 0.D
areac = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0    
    areac = !Pi * sma^2  - areac0 ;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*(intens-2.)
        fluxc0 = areac*(intens-2.)
        sma0 = sma
        ellip0 = ellip

    ENDIF

    IF area GT 0 THEN begin
        sum = sum + area*(intens-2.)
        sumc = sumc + area*(intens-2.)
        area0 = !Pi * sma^2*(1.-ellip)
        areac0 = !Pi * sma^2    ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0
;sumc = sumc - fluxc0
;mc = ((24.6 - 2.5*alog10(sumc)) - 39.26) 
;lc = 100^((4.63-mc)/5)

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclrtab2', /GET_LUN 

;read in the radial profile
rows= 18-2
rradius2 = FLTARR(rows)
rcounts2 = FLTARR(rows)
rerr2 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius2(j) = r
      rcounts2(j) = c
      rerr2(j) = ellerr
ENDFOR

close, lun
free_lun, lun

rcounts2 = rcounts2 - 2;[0:j-3] - 2
rradous2 = rradius2;[0:j-3]

result = MPFITFUN('exponential',rradius2[0:15],rcounts2[0:15], rerr2[0:15], start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


m = ((24.6 - 2.5*alog10(sum)) - 39.26) 
l = 100^((4.63-m)/5)

;----------------
openr, lun,"/n/Godiva1/jkrick/A3888/final/samemask/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*(intens-2.)
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*(intens-2.)

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclVtab2', /GET_LUN  ;iclVtab2

;read in the radial profile
rows=21
radius2 = FLTARR(rows)
counts2 = FLTARR(rows)
Verr2 = FLTARR(rows)
stop2 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
i = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      IF c GT 2.0 THEN BEGIN
          radius2(i) = r
          counts2(i) = c
          Verr2(i) = ellerr 
          stop2(i) = s
          i = i +1
      ENDIF
ENDFOR

close, lun
free_lun, lun

counts2 = counts2[0:i - 2] -  2
radius2 = radius2[0:i-2]

result = MPFITFUN('exponential',radius2[1:16],counts2[1:16], Verr2[1:16], start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

m = ((24.3 - 2.5*alog10(sum)) - 39.26) 
lB = 100^((4.82-m)/5)

;---
print, "A3888", lB , 100^((4.82-(((24.3 - 2.5*alog10(33922.8)) - 39.26)))/5),$
l, 100^((4.63-(((24.6 - 2.5*alog10(56704.4)) - 39.26)))/5)
;print, "A3888", l,lB,lc, l/lc
;-----------------------------------------------------------------

openr, lun,"/n/Godiva4/jkrick/A3984/elliptab", /get_lun

sum= 0
area0 = 0.D
areac0 = 0.D
areac = 0.D
area = 0.D
intens=0.D
i = 0
rvir = (34.9 *126.)

WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
    ellipf = ellip   ;final ellipticity
    smaf = sma       ;final sma
    areaf = area0    ;final ellipses area

ENDWHILE
close, lun
free_lun, lun

sum = sum - flux0

OPENR, lun,'/n/Godiva4/jkrick/A3984/iclrtab', /GET_LUN

;read in the radial profile
rows= 17  - 3
rradius = FLTARR(rows)
rcounts = FLTARR(rows)
rerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius(j) = r
      rcounts(j) = c 
      rerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun

result = MPFITFUN('exponential',rradius,rcounts, err, start, /quiet)   

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


;add to the sum from points outside of the end radius already summed
area0 = areaf
bins = 100.
FOR n = 1, bins , 1 DO BEGIN
    sma = (rvir- smaf)/bins * n + smaf
    area = !Pi * sma^2*(1.-ellipf) - areaf
    rad = sma - (rvir-smaf)/(2*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    areaf = area + areaf
;    print, "A3984",sma, area, intens, sum
ENDFOR


m = ((24.6 - 2.5*alog10(sum)) - 39.71) 
l = 100^((4.63-m)/5)

;-----------------
;r galaxies
openr, lun,"/n/Godiva4/jkrick/A3984/ellipgalrtab", /get_lun

sum= 0
area0 = 0.D
areac0 = 0.D
areac = 0.D
area = 0.D
intens=0.D
i = 0
rvir = (34.9 *126.)

WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
    areac0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
    ellipf = ellip   ;final ellipticity
    smaf = sma       ;final sma
    areaf = area0    ;final ellipses area

ENDWHILE
close, lun
free_lun, lun

sum = sum - flux0

OPENR, lun,'/n/Godiva4/jkrick/A3984/galrtab', /GET_LUN

;read in the radial profile
rows= 17  - 3
rradius = FLTARR(rows)
rcounts = FLTARR(rows)
rerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius(j) = r
      rcounts(j) = c 
      rerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun

result = MPFITFUN('exponential',rradius,rcounts, err, start, /quiet)   

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR


;add to the sum from points outside of the end radius already summed
area0 = areaf
bins = 100.
FOR n = 1, bins , 1 DO BEGIN
    sma = (rvir- smaf)/bins * n + smaf
    area = !Pi * sma^2*(1.-ellipf) - areaf
    rad = sma - (rvir-smaf)/(2*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    areaf = area + areaf
;    print, "A3984",sma, area, intens, sum
ENDFOR


mg = ((24.6 - 2.5*alog10(sum)) - 39.71) 
lg = 100^((4.63-mg)/5)

;------------------
close, lun
free_lun, lun
openr, lun,"/n/Godiva4/jkrick/A3984/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva4/jkrick/A3984/iclVtab', /GET_LUN

;read in the radial profile
rows= 18
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts = counts[0:j-1] ;- 2    
radius = radius[0:j-1]

result = MPFITFUN('exponential',radius,counts, Verr, start, /quiet)   

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

m = ((24.3 - 2.5*alog10(sum)) - 39.71) 
lB = 100^((4.82-m)/5)

;-------------
print, "A3984", lB , 100^((4.82-(((24.3 - 2.5*alog10(15236.)) - 39.71)))/5),lg, $
l, 100^((4.63-(((24.6 - 2.5*alog10(25328.)) - 39.71)))/5)

;-----------------------------------------------------------------
;A141 done manually with ellipses, flux.sxc from eapertures.prof.cl

VICL = 813.03;591.88
rICL = 2566.26;2002.6

Vgalaxies = 12862;10336
rgalaxies = 23122;18547

print, "A141",100^((4.82-(((24.3 - 2.5*alog10(VICL)) - 40.29)))/5),$
100^((4.82-(((24.3 - 2.5*alog10(Vgalaxies)) - 40.29)))/5),$
100^((4.63-(((24.6 - 2.5*alog10(rICL)) - 40.29)))/5),$
100^((4.63-(((24.6 - 2.5*alog10(rgalaxies)) - 40.29)))/5)

color = 24.3 - 2.5*alog10(VICL) - (24.6 - 2.5*alog10(rICL))
print, color

;----------------------------------------------------------------
openr, lun,"/n/Godiva7/jkrick/A114/original/elliptab", /get_lun

sum= 0
areac0 = 0.D
areac = 0.D
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
    areac = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
        areac0 = !Pi * sma^2    ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0 + 133.   ; from the clump

close, lun
free_lun, lun

OPENR, lun,'/n/Godiva7/jkrick/A114/original/icl2rtab', /GET_LUN

;read in the radial profile
rows= 18
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts = counts[0:j-1]    ;get rid of the points which have bad stop codes
radius = radius[0:j-1]

result2 = MPFITFUN('exponential',radius,counts, Verr, start, /quiet)    ;ICL

;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

m = ((24.6 - 2.5*alog10(sum)) - 41.04) 
l = 100^((4.63-m)/5)

;-----------------------------
openr, lun,"/n/Godiva7/jkrick/A114/original/ellipBtab", /get_lun

sum= 0
area0 = 0.D
area = 0.D
intens=0.D
i = 0
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, sma, intens, ellip
    area = !Pi * sma^2*(1.-ellip) - area0
;    area = !Pi * sma^2  - area0;test a circle
    IF i EQ 0 THEN BEGIN
        flux0 = area*intens
        sma0 = sma
        ellip0 = ellip
    ENDIF
    IF area GT 0 THEN begin
        sum = sum + area*intens

        area0 = !Pi * sma^2*(1.-ellip)
;        area0 = !Pi * sma^2   ;circle test
    ENDIF
    i = i + 1
ENDWHILE

sum = sum - flux0 + 170.  ;from the clump

close, lun
free_lun, lun


OPENR, lun,'/n/Godiva7/jkrick/A114/original/icl2Vtab', /GET_LUN

;read in the radial profile
rows= 10
Vradius = FLTARR(rows)
Vcounts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      Vradius(j) = r
      Vcounts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

Vcounts = Vcounts[0:j-1]    ;get rid of the points which have bad stop codes
Vradius = Vradius[0:j-1]

result2 = MPFITFUN('exponential',Vradius,Vcounts, Verr, start, /quiet)    ;ICL
;print, sma0, ellip0
;add to the sum from points inside of the cuttof radius already summed
area0 = 0
bins = 100
FOR n = 1, bins , 1 DO BEGIN
    sma = sma0/bins * n
    area = !Pi * sma^2*(1.-ellip0) - area0
    rad = sma - sma0/(2.*bins) ;at middle of bin
    intens = result(0)*exp(-rad/(result(1)))
    sum = sum + area*intens
    area0 = area + area0
;    print, sma, rad, area, intens, sum
ENDFOR

m = ((24.3 - 2.5*alog10(sum)) - 41.04) 
lB = 100^((4.82-m)/5)

;-----------
print, "A114", lB , 100^((4.82-(((24.3 - 2.5*alog10(3491.1)) - 41.04)))/5),$
l, 100^((4.63-(((24.6 - 2.5*alog10(6869.0)) - 41.04)))/5)
;-----------------------------------------------------------------
;A118 done manually with ellipses, flux.sxc from eapertures.prof.cl

VICL = 1005.23
rICL = 2722.8 

;VICL = 860.53
;ricl = 2360.75

Vgalaxies = 8034.;6633.6
rgalaxies = 17218.2;16541

print, "A118",100^((4.82-(((24.3 - 2.5*alog10(VICL)) - 41.04)))/5),$
100^((4.82-(((24.3 - 2.5*alog10(Vgalaxies)) - 41.04)))/5),$
100^((4.63-(((24.6 - 2.5*alog10(rICL)) - 41.04)))/5),$
100^((4.63-(((24.6 - 2.5*alog10(rgalaxies)) - 41.04)))/5)


color = 24.3 - 2.5*alog10(VICL) - (24.6 - 2.5*alog10(rICL))
print, color


END
