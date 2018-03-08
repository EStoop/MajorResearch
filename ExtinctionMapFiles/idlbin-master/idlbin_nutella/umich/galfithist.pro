PRO histogram
close, /all
colors = GetColor(/load, Start=1)
numoobjects = 200
subgals= replicate({object, centerx:0D, centery:0D, re:0D, sma:0D, smb:0D},numoobjects)

i = 0
junk = " "
re = 0.0
openr, lun, "/n/Godiva1/jkrick/A3888/galfit/re.log", /get_lun
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun,  junk
    c = strsplit(junk , " ", /extract) 
    re = float(c[3])
    xcen = float(c[1])
    ycen = float(c[2])
    subgals[i] = {object, xcen, ycen, re, 0, 0}
   i = i + 1

ENDWHILE
close, lun
free_lun, lun

openr, lun, "/n/Godiva1/jkrick/A3888/galfit/block.log", /get_lun
j = 0
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun,  x,y,a,b

    FOR j = 0, numoobjects - 1, 1 DO BEGIN
        IF (x LT subgals[j].centerx + 1. AND x GT subgals[j].centerx - 1. AND $
          y LT subgals[j].centery + 1. AND y GT subgals[j].centery -1.) THEN BEGIN
            subgals[j].sma = a
            subgals[j].smb = b
        ENDIF
    ENDFOR

ENDWHILE

close, lun
free_lun, lun


;re values below 0.5 are wrong anyway
FOR j = 0, numoobjects - 1, 1 DO BEGIN
    IF (subgals[j].re LE 0.5) THEN subgals[j]={object, 0., 0., 0., 0., 0.}
ENDFOR

;find average re and make a sigmaclip above that
av = mean(subgals.re)
FOR j = 0, numoobjects - 1, 1 DO BEGIN
    IF (subgals[j].re GT 3*sqrt(av)) THEN subgals[j]={object, 0., 0., 0., 0., 0.}
ENDFOR
rearr = fltarr(numoobjects)
smaarr = fltarr(numoobjects)
i = 0
FOR j = 0, numoobjects - 1, 1 DO BEGIN
    IF (subgals[j].re NE 0.0) THEN BEGIN
     rearr(i) = subgals[j].re
     smaarr(i) = subgals[j].sma
 ENDIF
i = i + 1
ENDFOR
rearr = rearr(0:i - 1)
smaarr = smaarr(0:i - 1)

print, rearr
err = dindgen(numoobjects) - dindgen(numoobjects) + 1         ;make a simple error array, all even errors
start = [.08]
result = MPFITFUN('linear',smaarr,rearr, err, start);,PARINFO=pi);, weights=weight)


!p.multi = [0, 0, 1]
ps_open, file = "/n/Godiva1/jkrick/A3888/galfit/revssma.ps", /portrait, xsize = 4, ysize = 4

plot, subgals.sma, subgals.re, psym = 2,title = "re vs. mask size",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "semi-major axis (pix)", $
ytitle = "re (pix)", yrange = [0,20]

print, result
print, result(0)

oplot, findgen(80), result(0)*findgen(80);, color = colors.blue
xyouts, 10, 18, 'best fit is that I mask out to 12re',  charthick =3, charsize = 0.9
xyouts, 10, 17, "but do I believe re?  NO"
;plot, histogram(reval, min = 1, max = 30),xrange = [1,30],xtitle = "re(pix)",$
;ytitle = "N", title = "Histogram of re values for a run through galfit w/o central galaxies",$
;thick = 3, charthick = 3, xthick = 3, ythick = 3

;plot, histogram(smaar), xtitle = "semi-major axis of the mask(pix)", ytitle = "N",$
;thick = 3, charthick = 3, xthick = 3, ythick = 3, title = "mask sizes for all central galaxies",$
;xrange = [1,30]

ps_close, /noprint, /noid
undefine, subgals
end
