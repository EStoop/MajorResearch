Pro getfinf, fieldstring, name, fsize
;
; (c) 1992 U.S. Government / National Aeronautics and Space Administration
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    GETFINF digests the output of getflds to get the field information 
;
;DESCRIPTION:
;	getfinf digests the output of getflds to get the field information
;	it returns the info in a SIZE()-like array
;
;CALLING SEQUENCE:
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param	I/O    type		description
;  fieldstrings I      string		output from getflds 
;  name		O      string		field name by itself
;  fsize	O      int array	like SIZE()
;					fsize(0) = dimensions, 0 if not array
;					fsize(1) = dimension 1...
;					fsize(n+1) = type code
;					fsize(n+2) = total number of elements
;					type codes: 0 undefined
;						    1 byte
;						    2 integer
;						    3 long
;						    4 float
;						    5 double  
;						    6 complex *
;						    7 string
;						    8 structure *
;						    * - not used here
;
;WARNINGS:
;	This software assumes that the strings are in the format that
;	getflds returns.
;
;	IT IS EXTREMELY UNFORTUNATE THAT THE READ_RDF ROUTINE 
;	DOES NOT PASS CODES RATHER THAN STRINGS.  CODES, LIKE THOSE
;	USED BY IDL'S SIZE FUNCTION, WOULD HAVE BEEN MUCH MORE USEFUL,
;	AND WOULD HAVE ELIMINATED THIS ROUTINE!
;	SEE ALSO GETFLDS.PRO, WHICH RETRIEVES THE STRINGS.
;
;EXAMPLE: 
;
;  example one
;
;	fieldlist = getflds( 'dataset' )
;	print, fieldlist(0:1)
;	FIELD0 - Floating Point  (1:5,1:6)
;	FIELD1 - Integer (1:5)
;	getfinf, fieldlist(0), name, fsize
;	print, name, fsize(0:3)
;	FIELD0	2	4	5	6
;	getfinf, fieldlist(1), name, fsize
;	print, name, fsize(0:3)
;	FIELD0	1	3	5 	0
; ____
;
;  example two
;
; ____
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;MODIFICATION HISTORY
;	Created 1-Jan-1993 Peter Kryszak-Servin SPR 10177
;
;.TITLE
;getfinf
;-
;


name = strmid( fieldstring, 0, strpos(fieldstring,' ') )
tend = strpos( fieldstring, ' ',40)

if tend le 0 then begin
	type = strmid( fieldstring, 40, strlen(fieldstring) )
	endif $
else begin
	type = strmid( fieldstring, 40, tend-40 )
	endelse
	
dim  = strmid( fieldstring, 64, strlen(fieldstring)-64 )
lpos = strpos( dim, '(' ) 
rpos = strpos( dim, ')' ) - 1
dim  = strmid( dim, lpos+1, rpos-lpos) 
	
ndims = 0
minelem = 0 ; minimum number of elements
if dim ne '' then begin
	dims = StrParse( dim, delim=',')
	ndims = n_elements( dims )
	minelem = 1
	;print, 'ndims -> ', ndims,dims
        endif


fsize = lonarr( ndims + 3 )
fsize( ndims + 2 ) = minelem
fsize(0) = ndims
if ndims gt 0 then begin
   for i = 1, ndims do begin
	srange = strparse( dims(i-1), delim=':' )
	;print, 'srange -> ', srange
	fsize(i) = long( strtrim(srange(1),2) )
	fsize( ndims + 2 ) = fsize( ndims + 2 ) * fsize(i)
	endfor
   endif


; determine the type, see CSDR$SOURCE:[UFC]Read_RDF.C for more info 
CASE 1 of
	(strpos( type, 'Byte'       ) ne -1) : fsize(ndims+1) = 1
	(strpos( type, 'Word'       ) ne -1) : fsize(ndims+1) = 2
	(strpos( type, 'Longword'   ) ne -1) : fsize(ndims+1) = 3
	(strpos( type, 'F_floating' ) ne -1) : fsize(ndims+1) = 4
	(strpos( type, 'D_floating' ) ne -1) : fsize(ndims+1) = 5
	(strpos( type, 'complex'    ) ne -1) : fsize(ndims+1) = 6
	(strpos( type, 'Character'  ) ne -1) : fsize(ndims+1) = 7
	else : message, /continue, 'Not a supported type.'
	endcase


;print, 'fsize -> ', fsize
;pause

RETURN
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


