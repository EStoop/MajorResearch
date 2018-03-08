function ellipprof, filename, number,xcenter, ycenter, data, blankdata, maxradius

print, "filename", filename
colors = GetColor(/Load, Start=1)

radius = fltarr(1000)
inten = fltarr(1000)
error = fltarr(1000)
rad = 0.D
i = 0.D
shit = string(0.0)

length = 0.

;must get rid of all "indef" in the tab file
OPENR, lun,filename, /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "file did not open"

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, rad, i, shit
    IF (shit EQ 'INDEF') THEN shit = 0.0
    radius(length) = rad
    inten(length) = i
    error(length) = float(shit)

    length = length +1
ENDWHILE

;shorten arrays
radius = radius(0:length-1)
inten = inten(0:length-1)
error = error(0:length-1)

;curve fitting
;------------------------------------------------------------------------------
maxfit = 32
IF (length GT 37 )THEN maxfit = 37

err = dindgen(length) - dindgen(length) + 1 
start = [0.2,7.0]
result = MPFITFUN('devauc',radius(21:maxfit),inten(21:maxfit), err(21:maxfit), start, bestnorm = chi2)

print, "deVauc result", result;, "for the galaxy", xcenter, ycenter
;plot as a check on the profile
;!p.multi = [number, 3, 3]
name = strarr(1)
name = string(xcenter+1) + STRING(ycenter+1)

x = findgen(30)
oplot, radius, inten, linestyle = 5;, title = name, xrange = [0,30], yrange = [0,8], xtitle = 'pixels', ytitle = 'counts/s'
oplot, x, (result(0)) * (exp(-7.67*(((x/(result(1)))^(1.0/4.0)) - 1.0))), thick = 3, color = colors.green
name = string(chi2)+ string(result(0)) + string (result(1))
xyouts, 0., 6., name, charsize = 0.5, color = colors.green

;make a model of the galaxy and add it to a blank image
;return the galaxy model in array newgal

newgal = subdev(xcenter, ycenter, maxradius, blankdata, result)

data = data - newgal

close, lun
free_lun, lun


end
