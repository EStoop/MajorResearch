;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;June 2002
;Jessica Krick
;
;  This program interpolates over bad pixles.  It is used specifically
;  in making a blank sky flat field.  To do the interpolation it makes
;  a histogram of the values in a 1300 pixel area box around the bad
;  pixel, and then selects a replacement pixel value at random from
;  the histogram.
;
;  input: file with image names to be interpolated on
;  output: new image file with bad pixels interpolated over
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro interpolate

close, /all
filename = "hello"
xmax = 2047
ymax = 3147
hist = fltarr(1225)
seed = 2320

fullarray = 1296
OPENR, lun1, '~/idlbin/interpolate.sep00.input', /GET_LUN
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, filename
    print, "filename", filename
;filename="/n/Godiva4/jkrick/aug98/bsflatr.bpm"

    infile = filename + '.fits'
    outfile = filename + '.int.fits'
    print, "infile", infile
    ;read in the fits file 
    FITS_READ, infile, data, header
    print, "read file"

    ;for all the pixels, go through and
    ;find the ones that are marked as bad = 0
    FOR x = 0, xmax - 1, 1 DO BEGIN
        FOR y = 0 , 2900 - 1 , 1 DO BEGIN
            IF (data[x,y] LT 0.1) THEN BEGIN ;500 - 600
             ;messy around bad column
             ;              IF (x LT 900 AND x GT 1000 ) THEN BEGIN  ;1025
             ;                  print, "ignoring bad colum", x, y
             ;              ENDIF ELSE BEGIN
;            print, "working on", x, y
                xcenter = x
                ycenter = y
                
                xsize1 = 100    ;radius of the extraction box in xdir
                ysize1 = 100    ;radius of the extraction box in ydir
                xsize2 = 100    ;radius of the extraction box in xdir
                ysize2 = 100    ;radius of the extraction box in ydir
;;                    print, xcenter, ycenter, xsize, ysize
                    ;first find out if going to go off the edge
                    ;and keep that from happening
                IF (xcenter + xsize2 GT xmax) THEN BEGIN
;                        yplus = FLOAT(fullarray / (4 * (xmax - xcenter)))
;                        ysize = FIX(yplus + 1)
                    xsize2 = xmax - xcenter
                ENDIF
                    
                IF (xcenter - xsize1 LT 1) THEN BEGIN
;                        yplus = FLOAT(fullarray / (4* abs(xcenter-xsize)))
;                        ysize = FIX(yplus + 1)
                    xsize1 = xcenter
;                        print, "inside x lt 1, new xsize, ysize", xsize, ysize
                ENDIF
                
                IF (ycenter + ysize2 GE ymax) THEN BEGIN
;                       xplus = FLOAT(fullarray / (4 * (ymax - ycenter)))
;                       xsize = FIX(xplus + 1)
                    ysize2 = ymax - ycenter - 1
;                       print, "inside GT ymax, xsize, ysize", xsize, ysize
                ENDIF
                
                IF (ycenter - ysize1 LT 1) THEN BEGIN
;                       xplus = FLOAT(fullarray / (4* abs(ycenter-ysize)))
;                       xsize = FIX(xplus + 1)
                    ysize1 = ycenter 
                ENDIF
;                    print, xcenter, ycenter, xsize, ysize
                     ;then fill the array for the histogram
                i = double(0)
                hist = fltarr(400000)
                
                
;                    print, xcenter, ycenter, xsize1, xsize2, ysize1, ysize2
                
                FOR x1 = xcenter - xsize1, xcenter + xsize2,1 DO BEGIN
                    FOR y1 = ycenter - ysize1, ycenter + ysize2, 1 DO BEGIN
                        hist[i] = data[x1,y1]
                        i = i +1
                    ENDFOR
                ENDFOR
                
                                ;truncate the array to be fullarray in size
                trunchist = fltarr(fullarray)
                FOR j = 0, fullarray-1, 1 DO BEGIN
                    trunchist[j] = hist[j]
                ENDFOR
                
                
                                ;sort the histogram low to high
                sorthist = SORT(trunchist)
                newhist = hist[SORT(trunchist)]
                
                                ;pick a random number, make sure it is an integer
                random = RANDOMN(seed, /UNIFORM)
                random = FIX(random * fullarray)
                
                                ;don't want to pick a bad number
                WHILE (newhist[random] LE 0 ) DO BEGIN
                    random = RANDOMN(seed,/UNIFORM)
                    random = FIX(random * fullarray)
                ENDWHILE
                
                                ;change out the bad pixel for a new good one
                data[xcenter,ycenter] = newhist[random]
;                ENDELSE
                
            ENDIF
        ENDFOR
    ENDFOR
    
                                ;write out the new image
    print, "writing image"
    fits_write, outfile, data, header
    
ENDWHILE 

close, lun1                ;closes the file with the list of filenames
free_lun, lun1
END
