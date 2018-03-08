PRO SETUP_PS,OLD_DEVICE,PLOTFILE=PLOTFILE
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    SETUP_PS configures for production of landscape PS plot files.
;
;DESCRIPTION:
;    SETUP_PS sets the device correctly to produce landscape PS
;    plot files.  In addition, it returns the current device type
;    extant on the system when the routine was called, enabling the
;    user to reset the plot device without having specific knowledge
;    of its nature after the desired PS plot is produced.
;
;CALLING SEQUENCE:  
;    SETUP_PS [, OLD_DEVICE] [,PLOTFILE = Plotfile]
;
;ARGUMENTS: (I = input, O = output, [] = optional)
;    OLD_DEVICE     [O]  str     Name of IDL plot device at time
;				 routine was called.
;    PLOTFILE       [I]  str     Name of plot file to be created.
;                                Defaults to 'GENERIC.PS_PLOT'.
;
;WARNINGS:
;    When user if finished creating the desired PS plot, he/she must
;    then close the PS file and (if desired) reset the plot device
;    using the commands
;             DEVICE,/CLOSE
;             SET_PLOT, OLD_DEVICE
;
;EXAMPLES:
;    1.  To create a plot of the vector YY using scale XX and write
;    it to the PS file 'XX_YY.PS', then reset the the plot device:
;	SETUP_PS, Old_Dev, Plotfile = 'XX_YY.PS'
;          Plot, xx, yy, xstyle = 1, ystyle = 1, title = 'EXAMPLE'
;       DEVICE,/CLOSE
;       SET_PLOT, Old_Dev
;
;    2.  To create a plot of the vector YY using scale XX and write
;    it to the PS file 'GENERIC.PS_PLOT', then reset the the plot 
;    device:
;	SETUP_PS, Old_Dev 
;          Plot, xx, yy, xstyle = 1, ystyle = 1, title = 'EXAMPLE'
;       DEVICE,/CLOSE
;       SET_PLOT, Old_Dev
;
;    3.  To create a plot of the vector YY using scale XX and write
;    it to the PS file 'XX_YY.PS', 
;	SETUP_PS, Plotfile = 'XX_YY.PS'
;          Plot, xx, yy, xstyle = 1, ystyle = 1, title = 'EXAMPLE'
;       DEVICE,/CLOSE
;
;    4.  To create a plot of the vector YY using scale XX and write
;    it to the PS file 'GENERIC.PS_PLOT', 
;	SETUP_PS, Old_Dev 
;          Plot, xx, yy, xstyle = 1, ystyle = 1, title = 'EXAMPLE'
;       DEVICE,/CLOSE
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Finds the system variable for the device name then resets the 
;     device to PS and opens an appropriately named file for the plot.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     NONE
;
;MODIFICATION HISTORY:
;     Alice Trenholme, General Sciences Corp. February, 1991
;     Modified by ART 3/17/92 to use keyword.
;.TITLE
;Routine SETUP_PS
;-
	old_device = !d.name

	set_plot,'PS'
	if n_elements(PLOTFILE) eq 0 then PLOTFILE='generic.ps_plot'
	device,filename = PLOTFILE,/landscape

end
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


