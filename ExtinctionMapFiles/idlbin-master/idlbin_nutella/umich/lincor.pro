
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;April 2002
;
;This program finds and corrects for nonlinearities in the ccd
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lincor

close, /all		;close all files = tidiness

expname = strarr(1)
filename = "hello";strarr(1)  	;declare filename to be a stringarray
exptime = 0.0

xmax = 2048
ymax = 3148


OPENR, lun1, '/n/Godiva2/jkrick/sep00/listlincor', /GET_LUN

WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, filename
    print, 'working on', filename
                                ;add suffixes for appropriate input and output filenames
    imagefile = filename ;'/n/Godiva4/jkrick/aug98/aug20/'+filename + '.fits'
;    outfile = filename + '.l.fits'
    
    FITS_READ, imagefile, data, header
    
;    FOR x = 1,xmax-1,1 DO BEGIN
;        FOR y = 1,ymax-1,1 DO BEGIN
;            IF (data(x,y) LT 12725.0) THEN BEGIN
;;                z= 1 - (.98687 + 1.0058E-6*data(x,y))
 ;               z= 1 - (.9384 + 2.32E-6*data(x,y))

;                data(x,y) = data(x,y) + z*data(x,y)
;            ENDIF
;        ENDFOR
;    ENDFOR
    
    z= (.926 + 2.6299E-6*mean(data))
;    z= (.98776 + 5.765E-7*mean(data))
;    print, data[100,100],z[100,100]
    data = data/z
;    print, data[100,100]

   
    FITS_WRITE, imagefile, data, header
ENDWHILE

close, lun1
free_lun, lun1

END


;I wish I knew why I needed this? normalized = 69.900
