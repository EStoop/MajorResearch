PRO deleteme

close,/all


openr, lun, "/n/Godiva6/jkrick/oct98/oct21/lintestr2",/get_lun

transarr = fltarr(30)
rarr = fltarr(30)
i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    readf, lun, e, c
    rarr[i] = e
    transarr[i] = c
    i = i + 1
ENDWHILE

transarr = transarr[0:i-1]
rarr = rarr[0:i-1]

close, lun
free_lun, lun
    yaxis = rarr/transarr

plot, transarr, yaxis,psym = 2

err = dindgen(i) - dindgen(i) + 1
start = [1,60]
result = MPFITFUN('linear',transarr, yaxis, err, start)
oplot, findgen(300)*100,result(1) + result(0)*findgen(300)*100,$
  thick = 3
normalized = result(1) + result(0)*28000
print, "normalized", normalized
FOR n = 0,i-1,1 DO BEGIN
    yaxis(n) = yaxis(n) / normalized
ENDFOR
plot, transarr , yaxis, psym = 2,thick = 3
start = [1,60]
result = MPFITFUN('linear',transarr, yaxis, err, start)


close, /all
END
