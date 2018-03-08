Pro pause, time=time, prompt=prompt
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PAUSE is used to pause for user to confirm reading a screen.
;
;DESCRIPTION:
;	PAUSE is used to pause for the user to confirm going ahead.
;	It prints the message give by the prompt keyword or a default.
;	It proceeds upon any keystroke, but will wait only as long
;	as the number of seconds given by the time keyword.
;
;CALLING SEQUENCE:
;	pause, time=30	;waits 30 seconds or until a key is pressed
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param	I/O    type		description
;    time	[I]    numeric		maximum time to wait, seconds
;    prompt	[I]    string		prompt string
;
;WARNINGS:
;
;EXAMPLE: 
;
;	print, 'Lots of information...'
;	pause, time=60, prompt='Press something to continue'
;
;	This example will print the prompt and then wait for any
;	keypress, but it will only wait for sixty seconds.
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
; Uses SysTime(), a standard IDL function
;
;MODIFICATION HISTORY
;	Date		SPR	Programmer / Notes
;	15-Dec-1992 	10178	Created by Pete Kryszak-Servin
;
;.TITLE
;PAUSE 
;-
;

If ( not (Keyword_Set( prompt )) ) Then Begin
	prompt = 'Press Any Key To Continue...'
	EndIf
Print, prompt

If ( keyword_set(time) ) Then Begin
	time1 = systime( 1 )
	While ( get_kbrd(0) eq '' ) Do Begin
		time2 = systime( 1 )
		If ((time2 - time1) GE time ) Then Return
		EndWhile
	Return
	EndIf $
Else Begin
	tmp = get_kbrd(1) 
	Return
	EndElse
	
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


