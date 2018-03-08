$ save_ver = f$verify()
$ if "''f$mode()'" .eqs. "BATCH" then set noverify
$!
$!******************************************************************************
$!
$!	This command procedure invokes COBE specific IDL.
$!
$! PLEASE READ FOLLOWING DIRECTIONS:
$!
$! Edit file such that replace everywhere #CGIS_IDL# with the top
$! level location of the CGIS software:
$!       e.g.   #CGIS_IDL# --> idl_dir:[lib.cgis]
$!
$! Note, take off final delimiter (']') in replacement for reprojection
$! lookup table
$!        e.g.  #CGIS_IDL#.qlut] --> idl_dir:[lib.cgis.qlut]
$!
$!***********************************************************************
$ Set NoOn
$!
$! Define IDL_PATH logical name.
$ PATH = "+#CGIS_IDL#,+IDL_DIR:[LIB]"
$!
$!Define IDL_STARTUP logical name, saving the user's definition:
$ If "''F$TrnLNm("IDL_STARTUP")'" .nes. ""
$ then
$       Define/User_Mode/NoLog USER_START "''F$TrnLNm("IDL_STARTUP")'" 
$ Else
$       Define/User_Mode/NoLog USER_START #CGIS_IDL#AuxStart.Pro
$ EndIf
$!
$!define directory for reprojection lookup tables
$ Define/User_Mode/Nolog CGIS_DATA #CGIS_IDL#.qlut]
$!
$ Define/User_Mode/NoLog IDL_STARTUP #CGIS_IDL#STARTUP.PRO
$!
$ Define/User_Mode/NoLog uimage_help #CGIS_IDL#uimage_help.dat
$!
$!define directory for default IDL save sets
$ Define/User_Mode/Nolog CGIS_CISS ADBDISK:[CGIS_CISS]
$!
$!define directory for default FITS files
$ Define/User_Mode/Nolog CGIS_FITS ADBDISK:[PDS_FITS]
$!
$ LOGICAL = F$TRNLNM("IDL_PATH")
$ IF LOGICAL .NES. ""
$ THEN
$ 	INDEX = F$TRNLNM("IDL_PATH",,,,,"MAX_INDEX")
$ 	IF (INDEX .GT. 0)
$ 	THEN
$		COUNTER = 1
$		LOOP:
$		NAME = F$TRNLNM("IDL_PATH",,COUNTER)
$		LOGICAL = LOGICAL + "," + NAME
$		COUNTER = COUNTER + 1
$		IF COUNTER .LE. INDEX THEN GOTO LOOP
$		LOGICAL = F$EDIT(F$EXTRACT(0,F$LENGTH(LOGICAL),LOGICAL),"TRIM")
$	ENDIF
$	DEFINE/USER_MODE/NOLOG IDL_PATH "''LOGICAL'","''PATH'"
$ ELSE
$	DEFINE/USER_MODE/NOLOG IDL_PATH "''PATH'"
$ ENDIF
$!
$!Determine the current graphics device.
$ IF F$TRNLNM("DECW$DISPLAY") .EQS. ""
$ THEN
$       IF F$TRNLNM("IDL_DEVICE") .EQS. ""
$       THEN
$	    IF F$GETDVI ("SYS$COMMAND", "TT_REGIS") .EQS. "TRUE"
$	    THEN
$               DEFINE/USER_MODE/NOLOG IDL_DEVICE REGIS
$	    ELSE
$               DEFINE/USER_MODE/NOLOG IDL_DEVICE TEK
$	    ENDIF
$       ENDIF
$ ELSE
$	DEFINE/USER_MODE/NOLOG IDL_DEVICE X
$ ENDIF
$!
$! Invoke IDL.
$ DEFINE/USER_MODE/NOLOG SYS$INPUT SYS$COMMAND
$ idl
$!
$ if "''f$mode()'" .eqs. "BATCH" then toss = f$verify(save_ver)
$! return to original IDL state
$ deassign idl_startup
$ @idl_dir:[bin]idl_setup
$!
$!DISCLAIMER:
$!
$!This software was written at the Cosmology Data Analysis Center in
$!support of the Cosmic Background Explorer (COBE) Project under NASA
$!contract number NAS5-30750.
$!
$!This software may be used, copied, modified or redistributed so long
$!as it is not sold and this disclaimer is distributed along with the
$!software.  If you modify the software please indicate your
$!modifications in a prominent place in the source code.  
$!
$!All routines are provided "as is" without any express or implied
$!warranties whatsoever.  All routines are distributed without guarantee
$!of support.  If errors are found in this code it is requested that you
$!contact us by sending email to the address below to report the errors
$!but we make no claims regarding timely fixes.  This software has been 
$!used for analysis of COBE data but has not been validated and has not 
$!been used to create validated data sets of any type.
$!
$!Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.
$!
