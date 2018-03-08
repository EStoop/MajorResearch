Function getflds, filename 
;
; (c) 1992 U.S. Government / National Aeronautics and Space Administration
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    GETFLDS calls and digests the output of Read_RDF for the specified file.
;
;DESCRIPTION:
;	GetFlds calls Read_RDF for the specified filename and digests the
;	output so that a user can read it easily.  The output may be used
;	in a menu selection of field names.
;
;CALLING SEQUENCE:
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param	I/O    type		description
;  filename	I      scalar string	filename or dataset name
;  getflds	O      array string	fields and data types, etc
;
;WARNINGS:
;	The data set must be in the DAFS and have an RDL in the 
;	CSDR$RDL directory.  That means that this is useful only
;	for COBEtrieve data sets.  
;	
;	This uses the COBEtrieve data server.
;
;	This software assumes that any file specification is VMS
;	format.
;
;	IT IS EXTREMELY UNFORTUNATE THAT THE READ_RDF ROUTINE 
;	DOES NOT PASS CODES RATHER THAN STRINGS.  CODES, LIKE THOSE
;	USED BY IDL'S SIZE FUNCTION, WOULD HAVE BEEN MUCH MORE USEFUL.
;	SEE ALSO GETFINF.PRO, WHICH DECODES THE STRINGS.
;
;EXAMPLE: 
;
;  example one
;
;	fieldlist = getflds( 'dataset' )
;	print, fieldlist
;	FIELD1 - Scalar Integer
;	FIELD2 - Floating Point Array (1:5,1:6)
; ____
;
;  example two
;
;	tfield = getflds( 'dataset' )
;	fcount = n_elements( tfield )
;
;	fieldlist = strarr( fcount + 4)
;	fieldlist( 0 ) = 'Select Field'   ;first
;	fieldlist( fcount + 1 ) = '    '
;	fieldlist( fcount + 2 ) = 'Help'
;	fieldlist( fcount + 3 ) = 'Exit' ;last options
;
;	status = umenu( fieldlist, title = 0, init=1, xpos=0)
;	
;	; status is the user's choice
; ____
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;	trim filename down to data set name
;	call Read_RDF
;	trim down output from Read_RDF
;	return the strings
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;MODIFICATION HISTORY
;	14-Jan-1993	10177	Kryszak	 created
;
;.TITLE
;GetFlds
;-
;
; trim down the file name to data set name
;
;print, 'filename = ', filename
colon 	= StrPos( filename, ':' )
bracket = StrPos( filename, ']' )
start   = 1 + Max( [colon, bracket])
;print, 'colon = ', colon
;print, 'bracket = ', bracket
;print, 'start = ', start

length  = StrLen( filename )
rdlname = StrMid( filename, start, length )
;print, 'length = ', length
;print, 'filename = ', rdlname

semicol = StrPos( rdlname, ';' )
if semicol ne -1 then rdlname = StrMid( rdlname, 0, semicol )
;print, 'semicol = ', semicol 
;print, 'rdlname = ', rdlname
period  = StrPos( rdlname, '.' )
if period ne -1 then rdlname = StrMid( rdlname, 0, period )
;print, 'period = ', period
;print, 'rdlname = ', rdlname

; 
; call Read_RDF
print, 'Retrieving field names now.'
text = Read_RDF( rdlname )

fcount = N_Elements( text )
;print, 'fcount = ', fcount
;print, 'full text = '
;for i = 0,fcount-1 do print, '|',text(i),'|'


; don't list structures
index = where ( (strpos(text, 'TYPE:Structure' ) eq -1), count )
if count ne 0 then text = text( index )

; trim text of any blanks
text = StrTrim( text, 2 )

; don't list empty lines (READ_RDF sometimes returns these)
index = where ( text ne '' ) 
if index(0) ne -1 then text = text( index )


fcount = N_Elements( text )


name   = StrPos( text, 'LEN:') - 2		; end of name
type   = StrPos( text, 'TYPE:' ) + 5		; start of Type
class  = StrPos( text, 'CLASS:') - type	- 2	; end of Type
dim    = StrPos( text, 'DIM:' )			; start of Dim
length = StrLen( text ) - dim			; end of Dim

fpos = 0	; field name field 
tpos = 40	; type field
dpos = 65	; dimension field

For i = 0, fcount-1 Do Begin

   newtext = String( Replicate( 32B, 80 ) )

   StrPut, newtext, StrTrim( StrMid( text(i), 0, name(i)), 2)+' ', fpos
   StrPut, newtext, StrTrim( StrMid( text(i), type(i), class(i)), 2), tpos

   if dim(i) ne -1 then begin
 	StrPut, newtext, StrTrim( StrMid(text(i), dim(i)+4, length(i)),2),dpos
	endif

   text(i) = newtext

   EndFor


; trim it up again
text = StrTrim( text, 2 )

RETURN, text
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


