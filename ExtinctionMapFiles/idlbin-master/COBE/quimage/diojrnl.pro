PRO diojrnl, input_format, input_type, instrume,  dsname, title, faceno, $
  dsfield, subscr, output_format, outname
;+
;  DIOJRNL - a component of the DATAIO subsystem.  If journaling is
;  enabled (that is, the variable JOURNAL_ON has a value of 1), then
;  this routine will record, within the UIMAGE journal file, the
;  activities of the DATAIO software.
;#
;  Written by Pete Kryszak-Servin.
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;------------------------------------------------------------------------
  COMMON journal, journal_on, luj
  IF(NOT defined(journal_on)) THEN journal_on = 0
  IF(NOT journal_on) THEN RETURN
  PRINTF, luj, 'DATAIO subsystem'
  IF(input_type eq 'U') THEN BEGIN
     PRINTF, luj, '  user defined data set selected'
     IF(input_format eq 'CSM') THEN BEGIN
       PRINTF, luj, '  input format:  COBEtrieve  (used BUILDMAP)'
       GOTO, endreport
     ENDIF
    IF(input_format eq 'CISS') THEN BEGIN
      PRINTF, luj, '  input format:  COBE IDL Save Set  (used UPUTDATA)'
      GOTO, endreport
    ENDIF
  ENDIF
;
;  Record info about what was inputted.
;  ------------------------------------
  CASE input_format OF
    'CSM': BEGIN
	PRINTF, luj, '  input format:  COBEtrieve  (used READ_SKYMAP)'
	PRINTF, luj, '  data set contains ' + instrume + ' data'
	PRINTF, luj, '  file name:  ' + dsname
	PRINTF, luj, '  field name:  ' + dsfield + '(' + subscr + ')'
	IF((faceno ge 0) and (faceno le 5)) $
          THEN PRINTF, luj, '  face:  ' + STRTRIM(faceno, 2) $
          ELSE PRINTF, luj, '  entire cube selected'
	END
    'CISS': BEGIN
	PRINTF, luj, '  input format:  COBE IDL Save Set  (used RESTORE)' 
	PRINTF, luj, '  data set contains ' + instrume + ' data'
	PRINTF, luj, '  file name:  ' + dsname
        END
    'FITS': BEGIN
	PRINTF, luj, '  input format:  FITS  (used READFITS)'
	PRINTF, luj, '  data set contains ' + instrume + ' data'
	PRINTF, luj, '  file name:  ' + dsname
        PRINTF, luj, '  field name: ' + dsfield + '(' + subscr + ')'
	IF((faceno ge 0) and (faceno le 5)) $
          THEN PRINTF, luj, '  face:  ' + STRTRIM(faceno, 2) $
          ELSE PRINTF, luj, '  entire cube selected'
	END
    'UIMAGE': BEGIN
	PRINTF, luj, '  input format:  UIMAGE data-object'
	PRINTF, luj, '  data set contains ' + instrume + ' data'
	PRINTF, luj, '  data-object name:  ' + title 
        END
    ELSE: BEGIN
	PRINTF, luj, '  input format:  undefined' 
        GOTO, endreport
        END
  ENDCASE
;
;  Record info about what was outputted.
;  -------------------------------------
  CASE output_format OF
    'CSM': BEGIN
	PRINTF, luj, '  output format:  COBEtrieve  (unallowed)'
	GOTO, endreport
        END
    'CISS': BEGIN
	PRINTF, luj, '  output format:  COBE IDL Save Set'
	PRINTF, luj, '  output file name:  ' + outname
        END
    'FITS': BEGIN
	PRINTF, luj, '  output format:  FITS'
	PRINTF, luj, '  output file name:  ' + outname    
        END
    'UIMAGE': BEGIN
	PRINTF, luj, '  output format:  UIMAGE internal format'
	PRINTF, luj, '  data-object title:  ' + title 
        END
    ELSE: BEGIN
        PRINTF, luj, '  output format:  unsupported'
        END
  ENDCASE
endreport:
;
;  Close off the report with a row of hyphens.
;  -------------------------------------------
  PRINTF, luj, '----------------------------------------' + $
               '--------------------------------------'
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


