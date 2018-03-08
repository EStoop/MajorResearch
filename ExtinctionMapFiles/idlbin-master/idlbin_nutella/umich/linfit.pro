pro linfit

device, true=24
device, decomposed=0

colors = GetColor(/Load, Start=1)
close, /all		;close all files = tidiness

;set up for plotting
mydevice = !D.NAME
!p.multi = [0, 0, 2]
SET_PLOT, 'ps'

;device, filename = '/n/Godiva3/jkrick/sep9940/linearity.ps', /portrait, $
;                BITS=8, scale_factor=0.9 , /color

ps_open, file = '/n/Godiva3/jkrick/sep9940/linearity.ps',/por,ysize = 8,/color
crapx = [0,0]
crapy = [0,0]
;!p.region=[0,0.5,1,1]
plot, crapx,crapy, psym = 2,thick = 3, $
	XRANGE = [0,28000],yrange = [0.40,1.05],xtitle ='counts',$
	ytitle = 'normalized counts/second', title = 'Sep 1999 40inch site3',$
        ystyle = 1, xstyle = 1, charthick = 3, xthick = 3, ythick = 3

FOR j = 0, 2, 1 DO BEGIN

    names = ["/n/Godiva3/jkrick/sep9940/sep02/lintestr1","/n/Godiva3/jkrick/sep9940/sep02/lintestr2","/n/Godiva3/jkrick/sep9940/sep02/lintestr3"]
    
    openr, linlun, names(j),/get_lun
    rows= 50
    exparr = FLTARR(rows)
    countarr = FLTARR(rows)
    i = 0
    e = 0.
    c = 0.

    WHILE (NOT EOF(linlun)) DO BEGIN
        READF, linlun, e, c
        
        exparr(i) = e
        countarr(i) = c
        i = i + 1
    ENDWHILE
    close, linlun
    free_lun, linlun
    
    exparr = exparr[0:i-1]
    countarr = countarr[0:i-1]
    
    sortindex = Sort(exparr)
    sortedexptime = exparr[sortindex]
    sortedcounts = countarr[sortindex]

    yaxis = sortedcounts/sortedexptime
    
  IF j EQ 0 THEN begin  
      sortedexptime = sortedexptime(0:i-3)
      sortedcounts = sortedcounts(0:i-3)
      i = i - 3
  ENDIF

    err = dindgen(i) - dindgen(i) + 1
    start = [1,60]
    result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

    normalized = result(1) + result(0)*28000
    print, "normalized", names(j), normalized
;    normalized = yaxis(i-1)
    FOR n = 0,i-1,1 DO BEGIN
        yaxis(n) = yaxis(n) / normalized
    ENDFOR
    print, sortedcounts, yaxis

    oplot, sortedcounts , yaxis, psym = j+4,thick = 3

    start = [1,60]
    result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

    
    oplot, findgen(300)*100,result(1) + result(0)*findgen(300)*100,$
      thick = 3,color = colors.blue
;    oplot, sortedcounts, 0.98086 + (0.00000125)*sortedcounts, $
;      thick = 3, color = colors.red, linestyle = 2

    newstring = string(result(1)) + string(result(0))
    xyouts, 1.0E4, 0.59+ (j/30.), newstring, color = colors.blue


ENDFOR


;newstring2 = string(0.98086) + string(1.25E-6)
;xyouts, 1.0E4, 0.91, newstring2, color = colors.red


;############################

FOR j = 0, 2, 1 DO BEGIN

    names = ["/n/Godiva3/jkrick/sep9940/sep03/lintestr", "/n/Godiva3/jkrick/sep9940/sep03/lintestB","/n/Godiva3/jkrick/sep9940/sep03/lintest.after" ]

    openr, linlun, names(j),/get_lun
    rows= 50
    exparr = FLTARR(rows)
    countarr = FLTARR(rows)
    i = 0
    e = 0.
    c = 0.

    WHILE (NOT EOF(linlun)) DO BEGIN
        READF, linlun, e, c
        
        exparr(i) = e
        countarr(i) = c
        i = i + 1
    ENDWHILE
    close, linlun
    free_lun, linlun
    
    exparr = exparr[0:i-1]
    countarr = countarr[0:i-1]
    
    sortindex = Sort(exparr)
    sortedexptime = exparr[sortindex]
    sortedcounts = countarr[sortindex]

    yaxis = sortedcounts/sortedexptime
    
    IF j EQ 0 THEN begin
        sortedexptime = sortedexptime(0:i-3)
        sortedcounts = sortedcounts(0:i-3)
        i = i - 3
    ENDIF

    err = dindgen(i) - dindgen(i) + 1
    start = [1,60]
    result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

    normalized = result(1) + result(0)*28000
    print, "normalized", normalized
;    normalized = yaxis(i-1)
    FOR n = 0,i-1,1 DO BEGIN
        yaxis(n) = yaxis(n) / normalized
    ENDFOR
    print, sortedcounts, yaxis

    oplot, sortedcounts , yaxis, psym = 2,thick = 3, color = colors.red

    start = [1,60]
    result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

    
    oplot, findgen(300)*100,result(1) + result(0)*findgen(300)*100,$
      thick = 3,color = colors.red
;    oplot, sortedcounts, 0.98086 + (0.00000125)*sortedcounts, $
;      thick = 3, color = colors.red, linestyle = 2

    newstring = string(result(1)) + string(result(0))
    xyouts, 1.0E4, 0.56 - (j/30.), newstring, color = colors.red

endfor

;----------------------------------------------------------------
;!p.region=[0,0,1.0,0.5]
    plot, crapx,crapy, psym = 2,thick = 3, $
	XRANGE = [0,28000],yrange = [0.35,1.05],xtitle ='counts',$
	ytitle = 'normalized counts/second', title = 'Oct 1998 40inch site3',$
        ystyle = 1, xstyle = 1, charthick = 3, xthick = 3, ythick = 3

FOR j = 0, 1, 1 DO BEGIN

    names = ["/n/Godiva6/jkrick/oct98/oct21/lintestr1","/n/Godiva6/jkrick/oct98/oct23/lintest" ]

    openr, linlun, names(j),/get_lun
    rows= 50
    exparr = FLTARR(rows)
    countarr = FLTARR(rows)
    i = 0
    e = 0.
    c = 0.

    WHILE (NOT EOF(linlun)) DO BEGIN
        READF, linlun, e, c
        
        exparr(i) = e
        countarr(i) = c
        i = i + 1
    ENDWHILE
    close, linlun
    free_lun, linlun
    
    exparr = exparr[0:i-1]
    countarr = countarr[0:i-1]
    
    sortindex = Sort(exparr)
    sortedexptime = exparr[sortindex]
    sortedcounts = countarr[sortindex]
    

    yaxis = sortedcounts/sortedexptime
    
    print, names(j), " sortedexptime, sortedcounts,yaxis", sortedexptime, sortedcounts,yaxis

    err = dindgen(i) - dindgen(i) + 1
    start = [1,60]
    result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

    normalized = result(1) + result(0)*28000
    print, "normalized", normalized
;    normalized = yaxis(i-1)
    FOR n = 0,i-1,1 DO BEGIN
        yaxis(n) = yaxis(n) / normalized
    ENDFOR
    oplot, sortedcounts , yaxis, psym = j+4,thick = 3
    
    print, names(j), sortedcounts, yaxis

    start = [1,60]
    result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

    
    oplot, findgen(300)*100,result(1) + result(0)*findgen(300)*100,$
      thick = 3,color = colors.blue
;    oplot, sortedcounts, 0.98086 + (0.00000125)*sortedcounts, $
;      thick = 3, color = colors.red, linestyle = 2

    newstring = string(result(1)) + string(result(0))
    xyouts, 1.0E4, 0.56 - (j/30.), newstring, color = colors.blue

endfor

openr, lun, "/n/Godiva6/jkrick/oct98/oct21/lintestr2",/get_lun

countarr= fltarr(30)
exparr = fltarr(30)
i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    readf, lun, e, c
    exparr[i] = e
    countarr[i] = c
    i = i + 1
ENDWHILE

countarr = countarr[0:i-1]
exparr = exparr[0:i-1]

close, lun
free_lun, lun
    yaxis = exparr/countarr

;oplot, countarr, yaxis,psym = 2, thick = 3

err = dindgen(i) - dindgen(i) + 1
start = [1,60]
result = MPFITFUN('linear',countarr, yaxis, err, start)
;oplot, findgen(300)*100,result(1) + result(0)*findgen(300)*100,$
;  thick = 3
normalized = result(1) + result(0)*28000
print, "normalized", normalized
FOR n = 0,i-1,1 DO BEGIN
    yaxis(n) = yaxis(n) / normalized
ENDFOR
oplot, countarr , yaxis, psym = 2,thick = 3
start = [1,60]
result = MPFITFUN('linear',countarr, yaxis, err, start)
oplot, findgen(300)*100,result(1) + result(0)*findgen(300)*100,$
  thick = 3



;ps_close
device, /close
set_plot, mydevice
END
