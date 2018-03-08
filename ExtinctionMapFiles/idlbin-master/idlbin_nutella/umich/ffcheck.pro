PRO ffcheck
close, /all
openw, outlun, "/n/Godiva1/jkrick/A3888/backlist",  /GET_LUN
fits_read, "/n/Godiva1/jkrick/A3888/centerV2000.block.fits", data, header
data = data - 2.
sum = 0.
av = fltarr(5000)

i = 0
j = 0

FOR x = 10, 300, 115 DO BEGIN;300
    FOR y = 600,1885, 115 DO BEGIN;1950
        printf, outlun, x,x+114,y,y+114
        FOR subx = 0,114 DO BEGIN
            FOR suby = 0, 114 DO BEGIN
                IF data[x+subx,y+suby] GT -10 THEN BEGIN
;                    print, data[x+subx,y+suby]
                    sum = sum + data[x+subx,y+suby]
                    i = i + 1
;                    print,"x, y,  sum, i", x+subx,y+suby,sum,i
                ENDIF
            ENDFOR
        ENDFOR
        IF i GT 200 THEN BEGIN
            av(j) = sum / i
;            print, av(j), i
            j = j + 1
        ENDIF
        i = 0
    ENDFOR
ENDFOR

i = 0
FOR x = 10, 1410, 115 DO BEGIN;300
    FOR y = 1550,1885, 115 DO BEGIN;1950
        printf, outlun, x,x+114,y,y+114
        FOR subx = 0,114 DO BEGIN
            FOR suby = 0, 114 DO BEGIN
                IF data[x+subx,y+suby] GT -10 THEN BEGIN
                    sum = sum + data[x+subx,y+suby]
                    i = i + 1
 ;                   print,"x, y,  sum, i", x+subx,y+suby,sum,i
                ENDIF
            ENDFOR
        ENDFOR
        IF i GT 200 THEN BEGIN
            av(j) = sum / i
;            print, av(j), i
            j = j + 1
        ENDIF
        i = 0
    ENDFOR
ENDFOR

i = 0
FOR x = 600, 1885, 115 DO BEGIN;300
    FOR y = 10,440, 115 DO BEGIN;1950
        printf, outlun, x,x+114,y,y+114
        FOR subx = 0,114 DO BEGIN
            FOR suby = 0, 114 DO BEGIN
                IF data[x+subx,y+suby] GT -1 THEN BEGIN
                    sum = sum + data[x+subx,y+suby]
                    i = i + 1
 ;                   print,"x, y,  sum, i", x+subx,y+suby,sum,i
                ENDIF
            ENDFOR
        ENDFOR
        IF i GT 200 THEN BEGIN
            av(j) = sum / i
;            print, av(j), i
        j = j + 1
        ENDIF
        i = 0
    ENDFOR
ENDFOR

i = 0
FOR x = 1500, 1885, 115 DO BEGIN;300
    FOR y = 10,1885, 115 DO BEGIN;1950        
        printf, outlun, x,x+114,y,y+114
        FOR subx = 0,114 DO BEGIN
            FOR suby = 0, 114 DO BEGIN
                IF data[x+subx,y+suby] GT -1 THEN BEGIN
                    sum = sum + data[x+subx,y+suby]
                    i = i + 1
                ENDIF
            ENDFOR
        ENDFOR
        IF i GT 200 THEN BEGIN
            av(j) = sum / i
            j = j + 1
        ENDIF
        i = 0
    ENDFOR
ENDFOR


av = av[0:j-1]
print, j
;print, av

device, true=24
device, decomposed=0

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/Load, Start=1)

device, filename = '/n/Godiva1/jkrick/A3888/ffcheck.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


plothist, av, xhist, yhist, bin = 0.1,thick = 3, xrange = [-1,1]
;plot, xhist, yhist;, xrange = [-0.3,0.3], charthick = 3, xthick = 3, thick = 3,$
ythick = 3;, xtitle = "Intensity (counts/s)", title = "V-band Pixel Intensity Histograms", ystyle = 5

device, /close
set_plot, mydevice
close, outlun
free_lun, outlun


END
;openw, outlun, "/n/Godiva7/jkrick/galfit/expback/listresults",  /GET_LUN, /append
;printf, outlun, xcenter+1 ,ycenter+1, result
;close, outlun
;free_lun, outlun
