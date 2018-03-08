# Shell script to start up COBE IDL (UIDL)
#
# define your local setup
setenv	IDL_DIR		/usr/local/rsi/idl_5

#
# Define the CGIS environoment
setenv	CGIS		/usr/local/rsi/idl_5/cgis
setenv	CSRC		/usr/local/rsi/idl_5/cgis
setenv	csrc		$CSRC
setenv uimage_help      $CSRC/uimage_help.dat
setenv	CDATA		$CGIS/data
setenv	CGIS_DATA	$CGIS/qlut/	# note final '/'
setenv	CGIS_DOC	$CGIS/userdoc
setenv	CGIS_FITS	$CDATA/fits
setenv	CGIS_CISS	$CDATA/ciss

#

#
# Now define IDL_PATH to cover it all. Don't forget that
# the IDL Astronomy User's Library must be included in your
# IDL path, too. You may add it here if you wish.
setenv	IDL_PATH	"+$IDL_DIR/lib":"+$CSRC"

# Determine the terminal time and set IDL_DEVICE accordingly.
# Currently we test for the presence of the DISPLAY environment
# variable.  If it exists then set display to X otherwise set it to
# TEK.
if ($?DISPLAY == 0) then
    setenv IDL_DEVICE	TEK
else
    setenv IDL_DEVICE	X
endif

#
# Aliases
#alias	idl		$IDL_DIR/bin/idl
alias	idltool		$IDL_DIR/bin/idltool
alias   uimage          'idl "uim.pro"'

# Define the IDL startup file for UIDL.  First save current value of 
# IDL_STARTUP.
if ($?IDL_STARTUP == 1) then
    setenv	USER_START	$IDL_STARTUP
    setenv	IDL_STARTUP	$CSRC/startup.pro
else
    setenv	USER_START	$CSRC/auxstart.pro
    setenv	IDL_STARTUP	$CSRC/startup.pro
endif

if ($IDL_STARTUP == $USER_START ) then
    setenv USER_START $CSRC/auxstart.pro
endif
#
# This start up widgets idl. If this is not defined on your system 
# make the appropriate adjustment.
idl

# If you want to keep this COBE IDL environment setup after executing
# this command comment out the next two lines:
# You may need to change the path depending on your local setup
unsetenv IDL_STARTUP
source $IDL_DIR/bin/idl_setup

#DISCLAIMER:
#
#This software was written at the Cosmology Data Analysis Center in
#support of the Cosmic Background Explorer (COBE) Project under NASA
#contract number NAS5-30750.
#
#This software may be used, copied, modified or redistributed so long
#as it is not sold and this disclaimer is distributed along with the
#software.  If you modify the software please indicate your
#modifications in a prominent place in the source code.  
#
#All routines are provided "as is" without any express or implied
#warranties whatsoever.  All routines are distributed without guarantee
#of support.  If errors are found in this code it is requested that you
#contact us by sending email to the address below to report the errors
#but we make no claims regarding timely fixes.  This software has been 
#used for analysis of COBE data but has not been validated and has not 
#been used to create validated data sets of any type.
#
#Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.
