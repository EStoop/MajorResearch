;+
; NAME:
;       DSNAME
; PURPOSE:
;       This function displays the input data set name menu.
; CALLING SEQUENCE:
;       DSNAME, input_format, dstype, dsname 
; INPUTS:
;       input_format     the input format, CSM, FITS, etc.
;       dstype           the input data set type, DIRBE, FIRAS, etc.
; OUTPUTS:
;  The data set name in text string format.
;       !Err = -1 if this function is called incorrectly
;	     = 88 is returned if 'Return to previous menu' is chosen
;
; COMMON BLOCKS:
;  None.
; RESTRICTIONS:
;  None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
; 17-Dec-1992	10330	Kryszak 	user interface mods
; 		10226	Kryszak 	FITS header display
;		10207	Kryszak         shows directory spec
; 23-Dec-1992	10387	Kryszak		help screens added
; 30-Dec-1992	10403	Kryszak		default directory fixed
;			Ewing		found CD for default dir
;                       Bazell          Changed prolog, updated 
;                                       comments.
; 5-Feb-1992    10463   Bazell          Changed calling sequence
;               10477                   of lstdsnam, changed the
;                                       order of some routines,
;                                       made some code more robust
; 11-Mar-1993   10672   Ewing           change backward control options
; 13-Apr-1993   10819   Ward            changes to handle passing around
;					the fits extension information
; 02-Jun-93     11003   Newmark         Change !version.os to !cgis_os. 
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
; SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
;  SPR 11981 08-Nov-1994  Read in skymap_info files.  J. Newmark 
;-
;
Pro dsname, input_format, dstype, dsname, fits_exten

!err = 1
!error = 1
status = 0


; The first three letters will be used for any subsystem specific files.
; for CSM and FITS format files.
; DIRBE  files are supposed to start with a B
; FIRAS  files are supposed to start with a F
; DMR    files are supposed to start with a D
; ADB    files are supposed to start with an X
; user defined, who knows?

Case dstype Of
	'B': subsys='B'
	'D': subsys='D'
	'F': subsys='F'
        'X': subsys='X'
	Else: subsys = input_format
        EndCase

; directory syntax handling
IF ( !cgis_os EQ 'vms' ) THEN BEGIN; it's vms
   pref = ''
   sep = ':'
   eod = ']'   ;for later
   index_filter = 'INDEX.MAP'
   archive_filter = 'ARCHIVE.'
ENDIF

IF (!cgis_os eq 'unix' ) then begin 
   pref = '$'
   sep = '/'
   eod = '/'
   index_filter = 'index.map'
   archive_filter = 'archive.'
ENDIF

IF (!cgis_os eq 'windows') then begin 
   pref = ''
   sep = '\'
   eod = '\'
   index_filter = 'index.map'
   archive_filter = 'archive.'
ENDIF


dirspec:
; The default directory spec is one of the following
;	 CSDR_ADB which includes the ADBs
;	 CGIS_FITS which points by default to CGIS$DATA
;	 CGIS_CISS which points by default to CGIS$DATA
; See UIDL.COM, its Unix equivalent, or the system environment
; variables for these values.  The startup procedures should
; allow the user to define these logical names / environment variables
;

Case input_format OF
	'CISS'	: default = 'CGIS_CISS'
	'FITS'	: default = 'CGIS_FITS'
	'CSM'	: default = 'CSDR_ADB'
	Else : goto, abort
	EndCase

; if the defaults are not defined then use current directory
lognam = getenv ( default ) 
If ( lognam Eq '' ) Then Begin
	cd, current=default
	lognam = default
	EndIf ; no log name

print, 'Enter the directory to search for your data.'
print, 'Press Return to use the default directory: '+default+'.'
print, 'Enter "current" to use your current directory.'
print, 'Enter "back" to return to the previous menu.'
If (default ne lognam) then print, default+' is defined as '+lognam+'.'
directory_spec = ''
read,  '> ', directory_spec

if STRUPCASE(STRTRIM( directory_spec, 2 )) eq 'BACK' then begin
	!ERR = -1
	!ERROR = -1
	goto,abort
	endif

if STRUPCASE(STRTRIM( directory_spec, 2 )) eq 'CURRENT' then $
    cd,current=directory_spec

; use the default if no directory spec entered
if (directory_spec eq '') then begin
	; if unix then prefix needed for environment variable
	; if vms then prefix is empty
	directory_spec = pref + default
        ; if windows then DOS doesn't know environment
        if (!cgis_os eq 'windows') then directory_spec=lognam
;        endif
  ENDIF ELSE lognam=directory_spec
;
; fix up the directory specification
; especially when the user entered something

directory_spec = STRTRIM( directory_spec, 2 )
; add directory separator if necessary
specend = STRPOS( directory_spec, eod)
specsep = STRPOS( directory_spec, sep)
if (specend eq -1) and (specsep eq -1) then begin
	directory_spec = directory_spec + sep
	endif
;
; check last character for proper path
;
IF (!cgis_os eq 'unix') OR (!cgis_os eq 'windows') THEN BEGIN 
  len1=strlen(lognam)
  last=strmid(lognam,len1-1,1)
  IF (last NE eod) THEN lognam=lognam+eod
  len1=strlen(directory_spec)
  last=strmid(directory_spec,len1-1,1)
  IF (last NE eod) THEN directory_spec=directory_spec+eod
ENDIF
;
; For CSM we assume that the data set starts with subsystem name.
; if user type of CSM, then full file spec required from buildmap
; (see dconvert.pro)

; Note: FITS and CISS files will be filtered just before display so that
; only valid FITS or CISS files are displayed as appropriate.
;     next line removed to remove instrument filtering for fits files
;     'FITS' : file_spec = directory_spec+subsys+'*.fit*'
CASE input_format of
     'CISS' : file_spec = directory_spec+'*.*'
     'FITS' : file_spec = directory_spec+'*.fit*'
     'CSM'  : file_spec = directory_spec+subsys+'*.*'
     Else   : goto, abort
EndCase
	
; get directory listing
status = lstdsnam( input_format, file_spec, flist, extension, fcount )

IF ( status LT 0 ) THEN Begin
   Message, 'No Data sets were found in the specifed directory', /continue 
   Goto, Abort ; return if errors
   EndIf

;
; filter out odd files, note: because with FITS files we only scan for
; *.FIT* files we don't need to filter files for them.  We do need to
; build the file and extension arrays however.
;
If input_format Eq 'FITS' Then Begin
	if (fcount gt 0) then begin
           tfiles = WHERE(strpos(strupcase(flist),'PIXPERM') eq -1)
           flist  = flist(tfiles)
           extension = extension(tfiles)
           tfiles = WHERE(strpos(strupcase(flist),'_MOD') eq -1)
	   tfile  = flist(tfiles)
           ext = extension(tfiles)
        endif else tfile = ['']
Endif

If input_format Eq 'CSM' Then Begin
	index = WHERE( StrPos( flist, index_filter ) Eq -1, cnt ) 
	if ( cnt gt 0 ) then tfile = flist( index ) 
EndIf

If input_format eq 'CISS' Then tfile = flist

fcount = n_elements( tfile )
if ( fcount eq 1 ) and ( strlen(tfile(0)) le 0 ) then fcount = 0

If fcount le 0 Then Begin
	print, 'No data sets were found in '+directory_spec
	print, ' '
	goto, dirspec
EndIf
lenl = strlen(lognam)
dirspecpos = 1 + STRPOS(tfile(0), eod, lenl-1 )
dirspec = STRMID( tfile(0), 0, dirspecpos )
file_list = strarr( fcount + 4)

file_list( 0 ) = 'Select Data Set - '+directory_spec ;first
file_list( fcount + 1 ) = '    '
file_list( fcount + 2 ) = 'HELP'
file_list( fcount + 3 ) = 'Return to previous menu'
file_list(1:fcount) = STRMID( tfile(0:fcount-1), dirspecpos, 255 )

menu:
	status = umenu( file_list, title = 0, init=1, xpos=0)
	If ( status LE 0 ) Then Begin 
	   !ERR = -1
	   !ERROR = -1
	   goto, abort ; return if errors
	EndIF

	if file_list(status) eq 'Return to previous menu' then goto, back

	if strupcase( file_list(status)) eq 'HELP' then begin
	    action, 'dsname.hlp' 
	    goto, menu
	EndIf

;	set data set name according to selection by user
	dsname = directory_spec + file_list( status )
	If (input_format Eq 'FITS') Then Begin
	    ; redo the extension list to start at 1, the same as
	    ; the file_list
	    file_ext  = strarr( fcount+1 )
	    file_ext(1:fcount) = ext
            nexten = n_elements(file_ext)
            if (status le nexten) then fits_exten = file_ext(status)
	    header = headfits( dsname )
	    review_header, header, 'FITS Header = '+dsname, use_file
	    if use_file eq 'N' then goto, back

	    ; skip over the binary header, as we don't have one in this case
	    index = where(strpos(header,'Written by UImage DATA IO') gt 0)
	    if index(0) ne -1 then goto, leaveok
            index = sxpar(header,'EXTEND')
            if (index ne 1) then goto, leaveok

	    ; read in the header for the fits BINARY TABLE extension
	    fxbopen,unit,dsname,1,bin_header
	    fxbclose, unit

	    ; setup use_file for review_header, then see if the user would
	    ; like to see it
	    use_file = 'N'
	    review_header, bin_header, 'FITS extension header', use_file, $
	           question='Would you like to review the Extension header?'

	EndIf 

leaveok:
	!Err = 1
	!Error = 1
	Return

back:
	!err = 88
	!error = 88
	dsname = ' '
	fits_exten = ' '
	Return

abort:
	Return

	End
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


