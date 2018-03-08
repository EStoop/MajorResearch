pro make_imask_from_dmask, dmaskfile, imaskfile, VERSION=version, $
                           BCDFILE=bcdfile
;+
; NAME:
;   MAKE_IMASK_FROM_DMASK
; PURPOSE:
;   Generate prototype imask from an input dmask
; EXPLANATION:
;   IMASK files are a cleaned version of pipeline produced dmask files; many
;   bits in the dmask have been set by the pipeline for internal purposes and
;   do not convey useful information to the observer.  This procedure creates
;   a new mask file by populating the desired bits from the dmask (crosstalk, 
;   radhit, latent, no flat applied, saturation and bad data bits) to the
;   new mask.  The new mask has placeholder bits for common IRAC data 
;   artifacts such as muxbleed, column pulldown and stray light while
;   remaining a 16 bit mask that can be used by MOPEX.  The imasks produced
;   have the same format and content as the imasks which will be produced by 
;   the SSC pipeline in future pipeline versions.
; NOTES:
;   Procedure requires astrolib procedures and functions to handle FITS I/O.
;   Output header could be cleaned up a little more.
;
; CALLING SEQUENCE:
;   MAKE_IMASK_FROM_DMASK, dmaskfile, [imaskfile, /VERSION]
;
; INPUT:
;   dmaskfile: string containing name of input dmask file
; OPTIONAL INPUT KEYWORD:
;   version: if set, then procedure prints version information and returns
;   bcdfile: string containing corresponding BCD filename, needed to add time
;             of observation to mask file
; OPTIONAL OUTPUT:
;   imaskfile: string containing name of output dmask file, if a filename is
;              not provided than the procedure will create one based on the
;              input filename.
;
; PROCEDURES CALLED:
;   WRITEFITS, SXADDPAR, SXADDHIST
; FUNCTION CALLED:
;   READFITS
;
; HISTORY:
;   Initial coding 19 Nov 2004 SJC
;   Made reading of FITS files mostly silent 22 Nov 2004 SJC
;   Added some header documentation and released 18 Mar 2005 SJC
;   Added option to input BCD file and write SCLK_OBS to imask file 26 April
;                                                      2005 SJC 
;   Corrected error in DBADBIT keyword 26 April 2005 SJC
;   Added handling of not linearized pixels 27 April 2005 SJC
;-

;	release_version = '1.0 18 MAR 2005'
	release_version = '1.1 26 APR 2005'
	mess = 'MAKE_IMASK_FROM_DMASK:' + !stime + ':'
	sclk_str = '[sec] SCLK time (since 1/1/1980) at DCE start'

	if (keyword_set(VERSION)) then begin
		print, 'MAKE_IMASK_FROM_DMASK Version ' + release_version
		return
	endif

	if (N_params() lt 1 or N_params() gt 2) then begin
		print, 'Syntax -- MAKE_IMASK_FROM_DMASK, dmaskfile[, imaskfile], $'
		print, '                                 BCDFILE=bcdfile'
		return
	endif

	if (N_params() eq 1) then make_imaskname = 1 else make_imaskname = 0

	dmask = readfits(dmaskfile, hdmask, /SILENT)
	himask = hdmask

; add info to header
	sxaddhist, mess + 'Creating imask', himask
	sxaddhist, 'from ' + dmaskfile, himask
; Describe bits in mask
	sxaddpar, himask, 'DQBIT', 0, ' Overall data quality represented by bit 0'
	sxaddpar, himask, 'ONEBIT', 1, ' Reserved for future use' 
	sxaddpar, himask, 'GHSTBIT', 2, ' Set if optical ghost present' 
	sxaddpar, himask, 'STRYBIT', 3, ' Set if stray light present' 
	sxaddpar, himask, 'DONTBIT', 4, ' Set if saturation donut' 
	sxaddpar, himask, 'MUXBIT', 5, ' Set if muxbleed present' 
	sxaddpar, himask, 'BANDBIT', 6, ' Set if banding present' 
	sxaddpar, himask, 'PULLBIT', 7, ' Set if column pulldown present' 
	sxaddpar, himask, 'XTLKBIT', 8, ' Set if crosstalk present' 
	sxaddpar, himask, 'RHITBIT', 9, ' Set if pixel contains radhit' 
	sxaddpar, himask, 'LATBIT', 10, ' Set if pixel contains latent image' 
	sxaddpar, himask, 'NOFLBIT', 11, ' Set if flat field was not applied latent image' 
	sxaddpar, himask, 'PLINBIT', 12, ' Set if pixel has large/no linearity correction' 
	sxaddpar, himask, 'SATBIT', 13, ' Set if pixel is saturated'
	sxaddpar, himask, 'DBADBIT', 14, ' Set if data is bad'


; Input bits to grab from dmask
	icrossbit = 2^13
	iradbit = 2^9
	ilatentbit = 2^5
	inotflatbit = 2^8
	isatbit = 2^10
	inonlinearbit = 2^12
	idbadbit = 2^14
	idmissbit = 2^11

; Bits to output for dmask set values
	ocrossbit = 2^8
	oradbit = 2^9
	olatentbit = 2^10
	onotflatbit = 2^11
	ononlinearbit = 2^12
	osatbit = 2^13
	odbogusbit = 2^14

; Initialize mask
	imask = dmask * 0

; Add crosstalk flagging
	ptr = where(dmask and icrossbit, count)
	if (count gt 0) then imask[ptr] = imask[ptr] or ocrossbit
	sxaddpar, himask, 'NCRSTLK', count, ' Number of pixels affected by crosstalk'

; Add radhit flagging
	ptr = where(dmask and iradbit, count)
	if (count gt 0) then imask[ptr] = imask[ptr] or oradbit
	sxaddpar, himask, 'NRADHIT', count, ' Number of radhit pixels'

; Add latent flagging
	ptr = where(dmask and ilatentbit, count)
	if (count gt 0) then imask[ptr] = imask[ptr] or olatentbit
	sxaddpar, himask, 'NLATENT', count, ' Number of pixels affected by latents'

; Add no flat field applied flagging
	ptr = where(dmask and inotflatbit, count)
	if (count gt 0) then imask[ptr] = imask[ptr] or onotflatbit
	sxaddpar, himask, 'NNOFLAT', count, ' Number of pixels with no flat-field correction'

; Add not linearized flagging
	ptr = where(dmask and inonlinearbit, count)
	if (count gt 0) then imask[ptr] = imask[ptr] or ononlinearbit
	sxaddpar, himask, 'NNONLIN', count, ' Number of pixels with no linearity correction'

; Add saturation flagging
	ptr = where(dmask and isatbit, count)
	if (count gt 0) then imask[ptr] = imask[ptr] or osatbit
	sxaddpar, himask, 'NSAT', count, ' Number of saturated pixels' 

; Add bad data flagging
	ptr = where((dmask and idbadbit) or (dmask and idmissbit), count)
	if (count gt 0) then imask[ptr] = imask[ptr] or odbogusbit
	sxaddpar, himask, 'NBAD', count, ' Number of bad pixels' 

; If output imask name is not provided, then make a reasonable one based
; on the input dmask name, that is swap the string "imask" for "dmask" or
; "bimsk" for "bdmsk"
	if (make_imaskname eq 1) then begin
		pieces = strsplit(dmaskfile, 'dmask', /EXTRACT, /REGEX)
		if (n_elements(pieces) eq 2) then $
		   imaskfile = pieces[0] + 'imask' + pieces[1] $
		else begin
			pieces = strsplit(dmaskfile, 'bdmsk', /EXTRACT, /REGEX)
			if (n_elements(pieces) eq 2) then $
			    imaskfile = pieces[0] + 'bimsk' + pieces[1] $
			else imaskfile = 'imask_' + dmaskfile
			
		endelse
	endif

	if (keyword_set(BCDFILE)) then begin
		himage = headfits(bcdfile)
		sclk = sxpar(himage, 'SCLK_OBS')
		sxaddpar, himask, 'SCLK_OBS', sclk, sclk_str
	endif

; Write imask
	writefits, imaskfile, imask, himask

return
end
