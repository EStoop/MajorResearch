;February 2005

;given a list of x y values for contours in and image
;(taken from ds9 contour)
;sum the counts inside of those contours


PRO readcontours
close,/all

histarr = fltarr(300000)
h = 0
fits_read, '/n/Godiva2/jkrick/A141/original/color.icl.fits', data, header
insidei = 0.D
sumi = 0.
insidej = 0.D
sumj = 0.
insidek = 0.D
sumk = 0.
insidel = 0.D
suml = 0.

FOR x = 1, 2848,1 DO BEGIN
    
    yarri = dblarr(4000000)
    i = 0
    columni = 0
    openr, luni, '/n/Godiva2/jkrick/A141/original/contouri', /get_lun
    WHILE (NOT eof(luni)) DO BEGIN
        readf, luni, xtest, ytest
        IF x LT xtest + 5. AND x GT xtest - 5. THEN BEGIN
            yarri[i] = ytest
            i = i + 1
            columni = 1
        ENDIF
        
    ENDWHILE
    ;----------------------------------------------------------------------
    yarrj = dblarr(4000000)
    j = 0
    columnj = 0
    openr, lunj, '/n/Godiva2/jkrick/A141/original/contourj', /get_lun
    WHILE (NOT eof(lunj)) DO BEGIN
        readf, lunj, xtest, ytest
        IF x LT xtest + 5. AND x GT xtest - 5. THEN BEGIN
            yarrj[j] = ytest
            j = j + 1
            columnj = 1
        ENDIF
        
    ENDWHILE
    ;----------------------------------------------------------------------
    yarrk = dblarr(4000000)
    k = 0
    columnk = 0
    openr, lunk, '/n/Godiva2/jkrick/A141/original/contourk', /get_lun
    WHILE (NOT eof(lunk)) DO BEGIN
        readf, lunk, xtest, ytest
        IF x LT xtest + 5. AND x GT xtest - 5. THEN BEGIN
            yarrk[k] = ytest
            k = k + 1
            columnk = 1
        ENDIF
        
    ENDWHILE
    ;----------------------------------------------------------------------
    yarrl = dblarr(4000000)
    l = 0
    columnl = 0
    openr, lunl, '/n/Godiva2/jkrick/A141/original/contourl', /get_lun
    WHILE (NOT eof(lunl)) DO BEGIN
        readf, lunl, xtest, ytest
        IF x LT xtest + 5. AND x GT xtest - 5. THEN BEGIN
            yarrl[l] = ytest
            l = l + 1
            columnl = 1
        ENDIF
        
    ENDWHILE
    ;----------------------------------------------------------------------

    IF columni EQ 1 OR columnj EQ 1 OR columnk EQ 1 OR columnl EQ 1 THEN begin
        yarri = yarri[0:i-1]
        sortindex= sort(yarri)
        sortedyarri = yarri[sortindex]
        num = N_ELEMENTS(sortedyarri) 
        
        FOR y = 1, 2878, 1 DO BEGIN
            IF y GE sortedyarri(0) AND  y LE sortedyarri(num-1) AND data[x,y] GT 0.0 $
                 AND data[x,y] LT 1.0 THEN BEGIN
                sumi = sumi + data[x,y]
                insidei = insidei + 1
            ENDIF            
        ENDFOR
        
    ENDIF

    IF (columnj EQ 1) AND (columnk EQ 0 ) AND (columnl EQ 0) THEN BEGIN
        yarrj = yarrj[0:j-1]
        sortindex= sort(yarrj)
        sortedyarrj = yarrj[sortindex]
        num = N_ELEMENTS(sortedyarrj) 
        
        FOR y = 1, 2878, 1 DO BEGIN
            IF y GE sortedyarrj(0) AND  y LE sortedyarrj(num-1) AND data[x,y] GT 0.0 $
                 AND data[x,y] LT 1.0 THEN BEGIN
                sumj = sumj + data[x,y]
                insidej = insidej + 1
            ENDIF            
        ENDFOR
        
    ENDIF
    
    IF columnk EQ 1 THEN begin
        yarrk = yarrk[0:k-1]
        sortindex= sort(yarrk)
        sortedyarrk = yarrk[sortindex]
        num = N_ELEMENTS(sortedyarrk) 
        
        FOR y = 1, 2878, 1 DO BEGIN
            IF y GE sortedyarrk(0) AND  y LE sortedyarrk(num-1) AND data[x,y] GT 0.0 $
                 AND data[x,y] LT 1.0 THEN BEGIN
                print, "inside k", x, y, data[x,y]
                histarr[h] = data[x,y]
                h = h + 1
                sumk = sumk + data[x,y]
                insidek = insidek + 1
            ENDIF            
        ENDFOR
        
    ENDIF
    
    IF columnl EQ 1 THEN begin
        yarrl = yarrl[0:l-1]
        sortindex= sort(yarrl)
        sortedyarrl = yarrl[sortindex]
        num = N_ELEMENTS(sortedyarrl) 
        
        FOR y = 1, 2878, 1 DO BEGIN
            IF y GE sortedyarrl(0) AND  y LE sortedyarrl(num-1) AND data[x,y] GT 0.0$
                 AND data[x,y] LT 1.0 THEN BEGIN
                suml = suml + data[x,y]
                insidel = insidel + 1
            ENDIF            
        ENDFOR
        
    ENDIF
    
    
    close, luni
    free_lun, luni
    close, lunj
    free_lun, lunj
    close, lunk
    free_lun, lunk
    close, lunl
    free_lun, lunl
    
    
ENDFOR

print, "sum", sumi, "inside"  , insidei , "mean  in i ", sumi / insidei
print, "sum", sumj, "inside"  , insidej, "mean  in j ", sumj / insidej
print, "sum", sumk, "inside"  , insidek , "mean  in k ", sumk / insidek
print, "sum", suml, "inside"  , insidel , "mean  in l ", suml / insidel

histarr = histarr[0:h-1]
plothist, histarr, thick = 3, bin = 0.01

END
