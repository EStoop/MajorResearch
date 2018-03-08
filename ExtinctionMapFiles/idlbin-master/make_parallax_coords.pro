pro make_parallax_coords, ra0, dec0, epoch0, pmra, pmdec, pi, epoch1, $
                          ra1, dec1, OBS0=obs0, OBS1=obs1, EDIR=edir, $
                          SILENT=silent, VERSION=version
                          
;+
; NAME:
;     MAKE_PARALLAX_COORDS
;
; PURPOSE:
;     Determine astrometric position in Equatorial (J2000) coodinates for an 
;     object given the position at a previous epoch, the proper motion and
;     the annual parallax.  Designed for refining coordinates for Spitzer 
;     observations.
;
;    
; CALLING SEQUENCE:
;     make_parallax_coords, ra0, dec0, epoch0, pmra, pmdec, pi, epoch1, $
;                           ra1, dec1, [OBS0=obs0, OBS1=obs1, EDIR=edir, $
;                           /SILENT]
;
; INPUTS:
;     ra0     - Right ascension of target for previous known epoch in degrees, 
;               scalar/vector float or double
;     dec0    - Declination of target for previous known epoch in degrees, 
;               scalar/vector float or double
;     epoch0  - Epoch of known coordinate in decimal years (e.g. 2004.5), 
;               scalar/vector float or double.  Acceptable values are between 
;               2003.66 and 2017.0 (the start of the Spitzer misson and the end 
;               of latest Spitzer ephemeris) for Spitzer as an observer and
;               between 1993.66 and 2017.0 for the Earth as an observer.
;     pmra    - Proper motion of object in Right ascension in arcseconds per
;               year, scalar/vector float or double
;     pmdec   - Proper motion of object in Declination in arcseconds per 
;               year, scalar/vector float or double
;     pi      - Magnitude of annual parallax in arcseconds / AU, scalar/vector 
;               float or double
;     epoch1  - Epoch of desired coordinate in decimal years, scalar/vector  
;               float or double.  Has same range as epoch0
;
;     Note ra0, dec0, epoch0, pmra, pmdec, pi all must have the same length and
;     either (ra0, dec0, epoch0, pmra, pmdec, pi) are vectors or epoch1 is a 
;     vector.  That is, the procedure will calculate a list of positions for
;     one output epoch and one or more input multiplets of input astrometric 
;     information or multiple output epochs and a single set of input 
;     astrometric information.  Multiple input observers cannot be mixed at this
;     time.
;
; OUTPUTS:
;     ra1     - Right ascension of target for desired epoch and observer in 
;               degrees, scalar/vector float or double
;     dec1    - Declination of target for desired epoch and observer in 
;               degrees, scalar/vector float or double
;
; OPTIONAL KEYWORD INPUTS:
;     OBS0    - Location of observer for input coordinates (ra0, dec0).  
;               Acceptable values are 'Spitzer', 'Earth' or 'Sun'.  If not 
;               provided, then input coordinate is assumed to be from Spitzer.
;     OBS1    - Location of observer for output coordinates (ra1, dec1). 
;               Acceptable values are 'Spitzer' or 'Earth', defaults to 
;               'Spitzer'
;
; OPTIONAL CONTROL KEYWORDS:
;     EDIR    - If set, specifies the path to the IDL save files containing
;               the required ephemerides.  Otherwise, defaults to the 
;               current working directory
;     SILENT  - If set, then no diagnostics are output
;     VERSION - If set, then print version of code (currently 1.2)
; 
; METHOD:
;     Starting with the known and hopefully correct coordinates of the first
;     provided epoch, the corrections to the RA and Dec relative to the Sun at 
;     each epoch for parallax of the observer are estimated using the formalism 
;     of Kirkpatrick et al. (2011, ApJS, 197, 19):
;            ra_parallax = pi * (R0 dot W^),
;            dec_parallax = -pi * (R0 dot N^),
;     where pi is the maximum parallax magnitude in arcseconds/AU, R0 is the
;     position vector of the observer relative to the Sun in Cartesian 
;     Equatorial coordinates viewed by the Sun, W^ is the local West unit vector  
;     for the direction of the source as viewed by the observer, N^ is the local 
;     North vector for the direction of the source as viewed by the observer.
;     That is, 
;            N^ = -cos(ra)sin(dec) x^ - sin(ra)sin(dec)y + cos(dec)z^, 
;            W^ = sin(ra) x^ - cos(ra) y^,
;     and the unit direction vector pointing along the source direction
;            R^ = cos(ra)cos(dec)x^ + sin(ra)cos(dec)y^ + sin(dec)z^.
;     ^ indicates unit vector, x^, y^, z^ are Cartesian Equatorial coordinates.
;     Note we ignore the difference between (RA, Dec) as observed from the Sun 
;     and the observer in this calculation.  The equations for parallax reduce
;     to those in the Explanatory Supplement of the Astronomical Almanac 
;     (Seidelmann 2005):
;            ra_parallax = X sin(ra) - Y cos(ra),
;            dec_parallax = pi * (X cos(ra) * sin(dec) + Y sin(ra) * sin(dec) 
;                                  - Z cos(dec)),
;     where X, Y and Z are the Cartesian position of the observer (in AU).  The
;     parallax is calculated based on observer position for each epoch resulting
;     in pairs of coordinate offsets (ra_parallax0, dec_parallax0) and
;     (ra_parallax1, dec_parallax1)
;    
;     In addition, the RA and Dec are updated using the provided proper 
;     motion offset between the two epochs:
;            dpmra = pmra * (epoch1 - epoch0)
;            dpmdec = pmdec * (epoch1 - epoch0)
;
;     The RA, Dec at the desired epoch are then
; 	        dec1 = dec0 + (dpmdec + dec_parallax1 - dec_parallax0) / 3600
;	        ra1 = ra0 + (dpmra  + ra_parallax1 - ra_parallax0) / (cos(dec)*3600)
;
; NOTES:
;     Ephemeris information is only provided for Spitzer and the geocenter 
;     of the Earth.  Observer information is relative to the barycenter.  Other
;     observer frames could be added by additional ephemerides.  Additional 
;     corrections are necessary for diurnal parallax for Earth-based observers.
;     Please send any questions or problems to the Spitzer help desk 
;     (help@spitzer.caltech.edu)
;
; EPHEMERIS FILE STRUCTURE:
;     Each observer ephemeris is an IDL save file containing the position of
;     the observer in Cartesian Equatorial coordinates (J2000) relative to the 
;     Sun as a function of time in years. The coordinates are on a 
;     spacing of two hours for Spitzer four hours for the Earth so there is 
;     minimal error in interpolating to the coordinate at the desired epoch.   
;     The file contains four vectors:
;          sepoch - dates of coordinates in decimal years
;          sx - x position (r cos(ra) cos(dec) of observer as seen by the Sun
;               (r is the distance to the Sun in AU, RA and Dec are the 
;                Equatorial coordinates of the observer as seen by the 
;                Sun)
;          sy - y position (r sin(ra) cos(dec) of observer as seen by the 
;                Sun
;          sz - z position (r sin(dec)) of observer as seen by the Sun
;     For some reason, HORIZONS cannot provide an ephemeris relative to the 
;     Barycenter.  The difference in parallax between using the Sun and the 
;     Barycenter should be small.
; 
;
; EXAMPLE:
;    One has three targets with earlier Earth-based observations with RA, DEC,
;    proper motion and annual parallax information, and wants to find out the
;    Spitzer RA and DEC for these three targets at 2013.234, using silent mode.
;    The IDL save files with coordinate information are in '/home/user1/idl/'.
;    IDL> make_parallax_coords,[345.251734,234.234454,69.234544],$
;    [-21.346721,34.4775,2.34556],[2004.753,1999.344,2001.344],$
;    [0.020,0.045,0.002],[-0.011,-0.045,0.223],[1.35,0.95,1.11],2013.234,$
;    ra1,dec1,OBS0='Earth',OBS1='Spitzer',EDIR='/home/user1/idl/',/SILENT  
;
; FUNCTIONS USED:
;     INTERPOL(), SIN(), COS()
;      
; REVISION HISTORY:
;     Code written, SJC 05 Dec 2012
;     Bug fixes, added an example, SL 07 Dec 2012
;     Update ephemerides with dates through 2019, added version keyword
;        SJC, 25 Oct 2016
;
;-

; Output informational message if calling syntax is incorrect
	syntax_message =['make_parallax_coords, ra0, dec0, epoch0, pmra, pmdec, $',$
	                  'pi, epoch1, ra1, dec1, [OBS0=obs0, OBS1=obs1, $', $
	                  'EDIR=edir, /SILENT']
	if (N_params() ne 9) then begin
		print, 'Please check syntax -- '
		for i = 0, n_elements(syntax_message)-1 do print, syntax_message[i]
		return
	endif
	
	if (keyword_set(VERSION) and not keyword_set(SILENT)) then $
		print, 'make_parallax_coords: version 1.2'
	
; Check inputs to make sure dimensionality is okay.  Code can either calculate
; final coordinates for a single input epoch and a vector of output epochs or
; final coordinates for a vector of input astrometry and a single output 
; epoch
	nra0 = n_elements(ra0)
	ndec0 = n_elements(dec0)
	nepoch0 = n_elements(epoch0)
	npmra = n_elements(pmra)
	npmdec = n_elements(pmdec)
	npi = n_elements(pi)
	
; Check to make sure input astrometry has same number of values for each input	
	diff = [nra0, ndec0, nepoch0, npmra, npmdec, npi] - nra0
	ptr = where(diff ne 0, dcount)

	if (dcount gt 0) then begin
		mstr = 'ra0, dec0, epoch0, pmra, pmdec, pi must have same number of '
		mstr = mstr + 'elements'
		message, mstr
	endif
	
	nepoch1 = n_elements(epoch1)
	if (nepoch1 gt 1) then begin
		if (nra0 gt 1) then $
		    message, 'Either input epoch or output epoch must be a scalar'
	endif 
	
; Copying input to temporary internal variables so that we can make sure that
; vectors of size one are treated as scalars
	if (nra0 eq 1) then begin
		ira0 = ra0[0] & idec0 = dec0[0] & iepoch0 = epoch0[0]
		ipmra = pmra[0] & ipmdec = pmdec[0] & ipi = pi[0]
	endif else begin
		ira0 = ra0 & idec0 = dec0 & iepoch0 = epoch0
		ipmra = pmra & ipmdec = pmdec & ipi = pi	
	endelse
				
; Control output, default is to output diagnostics
	if (keyword_set(SILENT)) then verbose = 0 else verbose = 1
	
; Set directory for ephemeris files if necessary
	if (keyword_set(EDIR)) then dir = edir else dir = './'
	if (file_test(dir, /DIRECTORY) eq 0) then message, 'EDIR not a directory!'
; Ephemeris file names
; Original ephemerides were relative to Sun and stopped in 2016
	spitzer_efile = 'spitzer_3vec_position_from_sun.sav'
	earth_efile = 'earth_3vec_position_from_sun.sav'
; Updated ephemerides are relative to Barycenter and go out to 2019
	spitzer_efile = 'spitzer_3vec_position_from_sun_2019.sav'
	earth_efile = 'earth_3vec_position_from_sun_2019.sav'

; Parse observer keywords and check for proper format
	m0string = 'OBS0 not Spitzer, Earth or Sun'
	if (keyword_set(OBS0)) then begin
		if (size(obs0, /TYPE) eq 7) then iobs0 = strlowcase(obs0) $
		else message, m0string
		if ((iobs0 NE 'spitzer') and (iobs0 NE 'earth') $
		         and iobs0 NE 'sun') then message, m0string
	endif else iobs0 = 'spitzer'
	
	m1string = 'OBS1 not Spitzer or Earth'
	if (keyword_set(OBS1)) then begin
		if (size(obs1, /TYPE) eq 7) then iobs1 = strlowcase(obs1) $
		else message, m1string
		if ((iobs1 NE 'spitzer') and (iobs1 NE 'earth')) then message, m1string
	endif else iobs1 = 'spitzer'	
	
; Determine observer position for first epoch, linearly interpolate between
; datum
	if (iobs0 EQ 'spitzer') then restore, '/Users/jkrick/idlbin/spitzer_3vec_position_from_sun_2019.sav' $;dir + spitzer_efile
	else restore, dir + earth_efile
	
; Check to see if desired epoch is within dates of ephemeris
        bptr = where(iepoch0 lt min(sepoch) or iepoch0 gt max(sepoch), bcount)
        print, 'min, max',  min(sepoch) , max(sepoch)
	if (bcount gt 0) then message, 'Epoch 0 = ' + strn(iepoch0[bptr]) + $
	                     ' outside of available ephemeris' 

; Calculate parallax even if Sun is reference
	if (iobs0 ne 'sun') then begin
		sx0 = interpol(sx, sepoch, iepoch0)
		sy0 = interpol(sy, sepoch, iepoch0)
		sz0 = interpol(sz, sepoch, iepoch0)	
	endif else begin
		sx0 = 0.0D
		sy0 = 0.0D
		sz0 = 0.0D
	endelse
; No reason to output zero parallax	
; first format for verbose output
	fstr0 = '("' + iobs0 + ' position for input epoch ", F7.2, " = ", A31)'

	if (verbose and iobs0 ne 'sun') then begin
		sr0 = sqrt(sx0*sx0 + sy0*sy0 + sz0*sz0)
		spra0 = atan(sy0, sx0) * 180.D / !DPI
		ptr = where(spra0 lt 0., count)
		if (count gt 0) then spra0[ptr] = spra0[ptr] + 360.D
		spdec0 = asin(sz0 / sr0) * 180.D / !DPI
		for i = 0, nra0-1 do print, FORMAT=fstr0, iepoch0[i], $
		                             adstring(spra0[i], spdec0[i], 4)
	endif

; If observer is not the same between epochs, then get the correct ephemeris
	if (iobs0 ne iobs1) then $
		if (iobs1 eq 'spitzer') then restore, dir + spitzer_efile $
		else restore, dir + earth_efile
		
; Check to see if desired epoch is within dates of ephemeris
	bptr = where(epoch1 lt min(sepoch) or epoch1 gt max(sepoch), bcount)
	if (bcount gt 0) then message, 'Epoch 1 = ' + strn(epoch1[bptr]) + $
	                     ' outside of available ephemeris'
	
; Determine observer position for second epoch, linearly interpolate between
; datum
; first format for verbose output
	fstr1 = '("' + iobs1 + ' position for output epoch ", F7.2, " = ", A31)'

	sx1 = interpol(sx, sepoch, epoch1)
	sy1 = interpol(sy, sepoch, epoch1)
	sz1 = interpol(sz, sepoch, epoch1)
	if (verbose) then begin	
		sr1 = sqrt(sx1*sx1 + sy1*sy1 + sz1*sz1)
		spra1 = atan(sy1, sx1) * 180.D / !DPI
		ptr = where(spra1 lt 0., count)
		if (count gt 0) then spra1[ptr] = spra1[ptr] + 360.D
		spdec1 = asin(sz1 / sr1) * 180.D / !DPI	
		nspra1 = n_elements(spra1)
		for i = 0, nspra1-1 do print, FORMAT=fstr1, epoch1[i], $
		                             adstring(spra1[i], spdec1[i], 4)
	endif
	
; Now calculate parallax for each epoch
; Need sin and cos of ra and dec of target first, we are using the first 
; epoch as the reference epoch.  This should be okay changes in coordinate
; are small from nominal epoch
	cra0 = cos(!DPI * ira0 / 180.D)
	sra0 = sin(!DPI * ira0 / 180.D)
	cdec0 = cos(!DPI * idec0 / 180.D)
	sdec0 = sin(!DPI * idec0 / 180.D)
	
	ra_parallax0 = ipi * (sx0 * sra0 - sy0 * cra0)
	dec_parallax0 = ipi * (sx0 * cra0 * sdec0 + sy0 * sra0 * sdec0 $
	                          - sz0 * cdec0)
	if (verbose) then begin
		for i = 0, n_elements(ra_parallax0)-1 do begin
			print, FORMAT='("Epoch 0 RA parallax (arcsec) = ", D9.4)', $
			        ra_parallax0[i]
			print, FORMAT='("Epoch 0 Dec parallax (arcsec) = ", D9.4)', $
					dec_parallax0[i]
		endfor
	endif

	ra_parallax1 = ipi * (sx1 * sra0 - sy1 * cra0)
	dec_parallax1 = ipi * (sx1 * cra0 * sdec0 + sy1 * sra0 * sdec0 $
	                          - sz1 * cdec0)
	if (verbose) then begin 
		for i = 0, n_elements(ra_parallax1)-1 do begin
			print, FORMAT='("Epoch 1 RA parallax (arcsec) = ", D9.4)', $
			        ra_parallax1[i]
			print, FORMAT='("Epoch 1 Dec parallax (arcsec) = ", D9.4)', $
					dec_parallax1[i]
		endfor
	endif
	
; and calculate proper motion difference,
	dpmra = ipmra * (epoch1 - iepoch0)
	dpmdec = ipmdec * (epoch1 - iepoch0)	
	if (verbose) then begin
		for i = 0, n_elements(dpmra)-1 do begin
			print, FORMAT='("RA PM component (arcsec) = ", D9.4)', dpmra[i] 
			print, FORMAT='("Dec PM component (arcsec) = ", D9.4)', dpmdec[i] 
		endfor
	endif	
	
; Now get new coordinates, remember proper motion and parallax corrections 
; are in arcseconds
	dec1 = dec0 + dpmdec / 3600.D + (dec_parallax1 - dec_parallax0) / 3600.D
	ra1 = ra0 + $
	    (dpmra / 3600.D + (ra_parallax1 - ra_parallax0) / 3600.D) / cdec0
	              
return
end
