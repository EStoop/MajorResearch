;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;November 2001
;Jessica Krick
;
;This program takes SExtractor output, (using the parameters
;	set in /n/Whalen1/jkrick/ICL/default.sex), reads in 
;	the appropriate fits image, and masks all the objects
;	in the image as found by SExtractor.
;
;	I make round masks for the saturated stars, so that the
;	diffraction spikes are also masked.
;	The grow parameter is just like that in iraf, I want to
;	make sure that the whole star is really masked
;
;This program was designed to mask out all the objects in a blank
;	sky image to get a flatter night sky flat
;
;input: in the same directory need to have a file named
;	lisblankskyr with the names of all the images to be masked
;
;output: an image named ccd****.m.fits with the objects set 
;	equal to 1E-7
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro remove_obj

close, /all

;declare variables
o = 0.0          	;object number
x = 0.0   		;xcenter
y = 0.0			;ycenter
e = 0.0   		;ellipticity
a = 0.0   		;semi-major axis
b = 0.0			;semi-minor axis
m = 0.0   		;magnitude (un-corrected for zpt, X, ...)
j = 0.0 		;junk
f = 0.0			;flux
fwhm = 0.0              
bkgd = 0.0
filename = strarr(1)  	;declare filename to be a stringarray


;read the filenames from a file:
OPENR, lun1, '~/idlbin/trash', /GET_LUN
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, filename
    print, 'working on', filename
    ;add suffixes for appropriate input and output filenames
    imagefile = filename + '.fits'
    datafile = filename + '.cat'
    newimagefile = filename + '.s2.fits'	

     ;read in the fits file and display original image

    FITS_READ, imagefile, data, header
    ;erase	
    ;newimage = congrid(data,507.5,510) ;makes smaller for display on screen
    ;tv, hist_equal(newimage)    ;displays


     ;open the SExtractor output datafile (ccd*.cat)
    OPENR, lun, datafile, /GET_LUN
    ;read the data file into the variables until hit EOF
    count = 0                   ;in case I want to know
    WHILE (NOT EOF(lun)) DO BEGIN
        READF, lun, o, x, y, a, b, e, f, m, j, j, fwhm, pa, bkgd
        
        IF (x EQ 1992 OR x EQ 1009 OR x EQ 1010 OR x EQ 1011 OR x EQ 1012 OR x EQ  1013 OR x EQ  1014 OR x EQ 1015 OR x EQ 1016 OR x EQ 1017 OR x EQ 1018 OR x EQ 1038 OR x EQ 1042 OR x GE 2045) THEN BEGIN
            ;these are the bad columns, don't want to grow them
        ENDIF ELSE BEGIN
            ;convert pa from degres into radians
            pa = pa*(!PI/180)
            ;figure out what the increment and grow should be 
            ;based on flux -> size of the object
            ;to deal with saturated stars
            IF (f GE 200000 ) THEN  BEGIN
                increment = (2*!PI)/1000
                b = a
                IF (a LT 15) THEN grow = 13 ELSE  grow = 7
            ENDIF
            IF (f GE 30000.0 AND f LT 200000) THEN BEGIN 
                increment = (2*!PI)/1000
                b = a		;make them circles the size of diff. spikes
                IF (a LT 15) THEN grow = 7 ELSE  grow = 4
            ENDIF
            IF (f LT 30000.0 AND f GT 8000.0) THEN BEGIN
                increment = (2*!PI)/50
                grow = 7
            ENDIF
            IF (f LE 8000.0) THEN BEGIN
                increment =  (2*!PI)/40
                grow = 7
            ENDIF
            
            IF (a GE 300) THEN grow = .5
            IF (m GE 80.0) THEN  BEGIN
                grow = 1        ;for the satellite
                a = 1
                b = 1
            ENDIF

            ;grow the axes
            ;grow = 7	;old method = grow them all by seven
            a = grow*a
            b = grow*b
            IF (a GE 60) THEN print, "a greater than 60",x,y,a

            IF (x GE 1000 AND x LE 1050) THEN print, x , y, a, "after growing"

            ;calculate which pixels need to be masked
            ;to do this change to a new coordinate system (ang)
            ; 	then work arond the ellipse, each time
            ;	incrementing angle by some small amount
            ;	smaller increments for dimmer (smaller) objs.
            FOR ang = 0.0, (2*!PI),increment DO BEGIN
                xedge = a*cos(pa)*cos(ang) - b*sin(pa)*sin(ang)
                yedge = a*sin(pa)*cos(ang) + b*cos(pa)*sin(ang)

                ; first get the center pixel
                IF (x GT 0 AND x LT 2030) THEN BEGIN
                    IF (y GT 0 AND y LT 2040) THEN BEGIN
                        data[x,y] = -100
                        ;make sure the other pixels are not off the edge
                        ;then mask the little suckers
                        IF (xedge+x GT 0 AND xedge+x LT 2030) THEN BEGIN
                            IF (yedge+y GT 0 AND yedge+y LT 2040) THEN BEGIN
                                ; have to do four cases to get everything 
                                ;between x and xegde and y and yedge
                                IF ((xedge+x) LT x) THEN BEGIN
                                    If (yedge+y LT y) THEN BEGIN
                                        data[xedge+x:x,yedge+y:y] = -100
                                    ENDIF ELSE BEGIN
                                        data[xedge+x:x,y:yedge+y] = -100
                                    ENDELSE
                                ENDIF ELSE BEGIN
                                    If (yedge+y LT y) THEN BEGIN
                                        data[x:xedge+x,yedge+y:y] = -100
                                    ENDIF ELSE BEGIN
                                        data[x:xedge+x,y:yedge+y] = -100
                                    ENDELSE
                                ENDELSE
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
                
            ENDFOR


            count = count + 1 
        ENDELSE  
    ENDWHILE

    ;close the SExtractor file, and free up the lun so it can be used again
    close, lun
    free_lun, lun

    ;for ccd3126 bleedout columns
    ;data[1052:1073,1022:1522] = -100
    ;data[1052:1073,496:996] = -100
    
    ;display the new image with the masked pixels
    ;twoimage = congrid(data,507.5,510) ;makes smaller for display on screen
    ;tv, twoimage                ;displays

    ;write a new fits file with the masked pixels
    fits_write, newimagefile, data, header

ENDWHILE

close, lun1                     ;closes the file with the list of filenames
free_lun, lun1

end
