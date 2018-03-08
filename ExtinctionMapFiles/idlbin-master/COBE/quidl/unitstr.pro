Function unitstr, units_code 
;+
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;	UnitStr returns a string with units corresponding to the units code.
;
; DESCRIPTION: 
;       This function returns the string corresponding to the units code.  
;	These codes are used in the ADB header files.  
;
;	The returned string will be empty if the units code is out of range.
;	Units code 0 corresponds to mixed units ('mixed').
;	Units code 1 corresponds to no units ('none').
;
;	See the file CSDR$SOURCE:[XTR]XTR_UNITS.TXT for the definitive list.
; 
;	
; CALLING SEQUENCE:
;	units_string = UnitStr( units_code )
;
; ARGUMENTS: 
;       units_code	Input 	int 	code for units, see XTR_UNITS.TXT
;
; RETURNS:
;       units_string	Output  string  units
;					empty string if code is out of range
;
; WARNINGS:
;	!err / !error 	are not set even if the units code is out of range
;	check the string, if it is empty then the units code is out of range
;
; EXAMPLE:
;	UIDL> print, UnitStr( 51 )
;	W/m**2/Hz/sr
;	UIDL>
;
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;	None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, September 1992.
;			Credit to Anne Raugh and Jeff Mullins of HSTX for
;			the FORTRAN version of the strings and codes from
;			which this is taken.
;			SER 9906, SAT CCR 47
;
;-
;
units_string = [        'mixed',  $
			'none',   $ 
			'h',      $ 
			'min',    $ 
			's',      $ 
			'deg',    $ 
			'arcmin', $ 
			'arcsec', $ 
			'rad',    $ 
			'mag',    $ 
			'Hz',     $ 
			'mm',     $ 
			'um',     $ 
			'Angstrom', $ 
			'nm',     $ 
			'%',      $ 
			'Jy',     $ 
			'MJy/sr', $ 
			'K',      $ 
			'W/m**2/sr', $ 
			'K km/s', $
			'GHz',    $ 
			'MHz',    $ 
			'sr',     $ 
			'mK',     $ 
			'G',      $ 
			'T',      $ 
			'Wb',     $ 
			'm',      $ 
			'cm',     $ 
			'km',     $ 
			'au',     $ 
			'pc',     $ 
			'kpc',    $ 
			'Mpc',    $ 
			'km/s',   $ 
			'km/s/kpc', $ 
			'm/s',    $ 
			'erg',    $ 
			'erg/s',  $ 
			'J',      $ 
			'J/m**3', $ 
			'eV',     $ 
			'keV',    $ 
			'W',      $ 
			'cd',     $ 
			'cd/m**2', $ 
			'lm',     $ 
			'lux',    $ 
			'W/m**2', $ 
			'W/m**2/Hz', $ 
			'W/m**2/Hz/sr', $      
			'g',      $ 
			'kg',     $ 
			'Msun'    $ 
			]


numunits = size( units_string )

; check the units code for bounds
if (units_code lt 0) or (units_code ge numunits(1)) then begin
        !err = -1
	!error = -1
	return, ''
        endif

; return the string after trimming blanks
return, units_string( units_code )

END
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.


