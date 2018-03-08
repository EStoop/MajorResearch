pro irac_p_correct_pulldown_standalone_roberta_ch1, ifile, sfile, mfile, ofile, osfile, $
                             omfile, pthres, ft, backg, PTHRESH=pthresh, PLOG=plog, $
                             PULLNOISE_FACTOR=pullnoise_factor, $
                             VERBOSE=verbose, USESIG_WEIGHT=usesig_weight, $
                             NORDER=norder
;+
; NAME:
;   IRAC_P_CORRECT_PULLDOWN_STANDALONE
;
; PURPOSE:
;   Correct for instances of pulldown in a BCD
;
; INPUTS:
;   ifile : string containing BCD filename
;   sfile : string containing BCD uncertainty image filename
;   mfile : string containing BCD imask filename
;
; OUTPUTS:
;   ofile : string containing output (corrected) BCD filename
;   osfile : string containing updated BCD uncertainty image filename
;   omfile : string containing properly pulldown masked BCD imask filename
;   If these keywords are not provided, then the input files are overwritten
;
; CONTROL KEYWORDS:
;    PTHRES : Floating point, if set this is the threshold in MJy/sr * second 
;             for pulldown to occur.
;    PULLNOISE_FACTOR : Floating point, if set this is the value to inflate the
;                       uncertainty image by after pulldown correction
;    USESIG_WEIGHT : Boolean, if set, then use the sigma image for weighting 
;                    the data in polynomial fitting the pulldown
;    PLOG : String, if set output diagnostic information to the logfile 
;           specified
;    VERBOSE : Boolean, if set, output diagnostic images of the various 
;              correction steps as well as other information
;
; ALGORITHM:
;    BCD, uncertainty and mask images are input
;    The mask image has column pulldown flagging undone
;    All pixels above the specified threshold in the image are identified
;         The x and y coordinates of the above threshold pixels are identified
;         The number of unique columns with a trigger are identified
;         The mask is updated with column pulldown flagging for those columns
;    An estimate of the image without artifacts is made
;    A difference image image - estimate is made
;    For each column containing a trigger
;         The peak pixel of the column is identified
;         The column is segmented into above the peak and below the peak
;              For each segment:
;                 The good pixels (not flagged as artifact on NaNed) are 
;                     identified; other pixels are NaNed.
;                  A polynomial is fit to the good pixels in the difference 
;                     vector for each segment using MPFIT.
;                  The segment fit is placed in the correction vector
;          The correction for the trigger is the average of the adjacent pixels 
;               in the correction vector
;          The correction vector is applied to the image
;          The uncertainty image should be updated for this column (maybe not 
;               yet implemented).
;    The corrected image, uncertainty and mask are written to the output 
;       filenames.
;
;
; NOTES:
;   BCD file is overwrittent with pulldown corrected BCD, uncertainty image 
;   is updated with revised uncertainty and imask is updated with new column 
;   pulldown flagging unless output filenames are specified.  Need to write 
;   header documentation.
;
; PROCEDURES USED:
;   SXADDPAR, WRITEFITS, IRAC_P_MAKE_GWEIGHTS, IRAC_P_ESTIMATE_SKY
;
; FUNCTIONS USED:
;   READFITS, SXPAR, BSORT, UNIQ, MEDIAN, MOMENT, ROBUST_SIGMA, POLY_FIT, MPFIT, 
;   POLY1D_FUNC_FOR_MPFIT, POLY
;
; HISTORY:
;    06 April 2005 SJC Segemented from fit_artifacts procedure
;    02 May 05 SJC added noise updating
;    04 May 05 SJC Use ROBUST_SIGMA to determine the uncertainty in corrections 
;    10 May 2005 SJC corrected sigma updating, handled NaNs in poly_fit
;    12 May 2005 SJC corrected handling of NaNs when optional polynomial 
;                    fitting is used
;    29 Jun 2005 SJC corrected estimation of corrected sigma on source
;    06 July 2005 SJC added check in poly fit to see if there is enough data 
;    12 Aug 2005 SJC changed procedure name
;    22 August 2005 SJC Added more diagnostics to log files
;    12 March 2006 SJC Added NAN keyword to max function
;    10 March 2010 SJC modified module to act as standalone, incorporated 
;                   estimate_truth and pulldown flagging in module, updated to 
;                   use MPFIT, code cleaned and made more efficient
;    11 March 2010 SJC confirmed that code compiles and works pretty much, 
;                  tinkering with the fitting is still required.
;-


; Define constants
	name = 'IRAC_P_CORRECT_PULLDOWN_STANDALONE_ROBERTA_CH1' 
	syntax =['irac_p_correct_pulldown_standalone_roberta, image, uncertainty, mask,$', $
	         '                                         oimage, ounc, omask, pull_th, int_time, backgr']
	
; SCLK time of warm mission start
	warm_start_sclk = 927192851.D
	
; Threshold for muxbleed to start in MJy * s / sr
        pulldown_thresh = [600, 600]   ; corr7
; Bit set in imask for muxbleed
	pulldownflag = 2^7

; Width for Gaussian smoothing kernel (sigma = FWHM / 1.665) in pixels
	sgauss = [1.44, 1.43] / 1.665 / 1.22
; Bits set for bad pixels
	bad_pix_bits = long(total(2L^[3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]))
; Set bits for data not to have solution fit to, assume that straylight, 
; muxbleed, banding are alright for this case (so 3, 5, 6 are okay)
	bad_notpull_bits = long(total(2L^[8, 9, 10, 11, 12, 13, 14, 15]))

; set maximum number of iterations
	max_iter = 1
 	
; Test inputs	
	if (N_params() ne 3 and N_params() ne 9) then begin
		print, syntax
		return
	endif

; If only input filenames passed then use these files as output as well	
	if (N_params() eq 3) then begin
		ofile = ifile
		osfile = sfile
		omfile = mfile
	endif
	
; Should also test for presence of input files
     if (n_elements(ifile) lt 1 or n_elements(sfile) lt 1 or n_elements(mfile) lt 1) then message,'Not valid input files'

; Set polynomial order, default to second order
        if (not keyword_set(NORDER)) then norder = 2	
	
; Set scaling of noise for pixels affected by muxbleed
	if (not keyword_set(PULLNOISE_FACTOR)) then pullnoise_factor = 1.25

; Set flags for outputting diagnostics
	if (not keyword_set(VERBOSE)) then begin
		silent = 1 
		write = 0 
	endif else begin
		silent = 0
		write = 1
	endelse

for iter =0, max_iter -1 do begin

 print,'iter =', iter
 if (iter eq 0) then begin


; Read in image file
	image = readfits(ifile, h, SILENT=silent)
	if (n_elements(image) eq 1) then message, 'Not a valid BCD file'
; Should test to make sure it is a BCD, otherwise this procedure is not useful
	naxis = sxpar(h, 'NAXIS')
	if (naxis eq 3) then message, 'Code does not work for subarray BCDs'

; Obtain information about image from header
	ch = sxpar(h, 'CHNLNUM') & print,'ch =',ch
	nx = sxpar(h, 'NAXIS1')
	ny = sxpar(h, 'NAXIS2')
	nxy = nx * ny

; There is no  pulldown for channels 3 and 4, just exit for those channels
	if (ch ge 3) then message, 'Module does not work for channels 3+4', $
	                           NOWRITE=nowrite
	
; Read in sigma file and mask
	sigma = readfits(sfile, hs, SILENT=silent)
	if (n_elements(sigma) eq 1) then message, 'Not a valid uncertainty image'
	mask = readfits(mfile, hm, SILENT=silent)
	if (n_elements(mask) eq 1) then message, 'Not a valid mask file'
	
; Determine if data are warm mission or cold
	t = sxpar(h, 'SCLK_OBS')
	if (t gt warm_start_sclk) then warm_flag = 1 else warm_flag = 0
	
; Print warning that procedure is only for warm data 
	if (warm_flag lt 1) then $
	    message, 'No fit for cryo IRAC data -- developmental pulldown version'
	
; Determine frametime of image
	ft = sxpar(h, 'FRAMTIME')
	
; Determine threshold (MJy/sr) for pulldown trigger if not specified
        if (not keyword_set(PTHRESH)) then pthres = pulldown_thresh[ch-1] / ft

 endif else begin

        image = image
        sigma = sigma
        mask = mask

 endelse

	
; Determine threshold (MJy/sr) for pulldown trigger

        ptr_neg = where(image lt 0, ocount_n)

        if (keyword_set(PTHRESH)) then begin

         if (ocount_n ne 0) then begin

           xneg = ptr_neg mod nx
           yneg = ptr_neg / ny
           plothist,xneg, xhist, yhist,/noplot
           med_neg = median(yhist)
           mean_neg = (moment(yhist))(0)
           if (med_neg eq 0) then med_neg = mean_neg
           sig_neg = sqrt((moment(yhist))(1))
           if (med_neg lt  (nx*0.75)) then begin

            ptr_p = where(yhist ge 3.*med_neg, ocount_p) 

           endif else begin 
  
            ptr_p = where(yhist lt (med_neg-sig_neg), ocount_p)
 
           endelse

           if (ocount_p ne 0) then begin 

            xp = xhist(ptr_p)
            print,'xp =',xp
            sptr_p = bsort(xp)
            xp = xp[sptr_p]
            sptr_p = uniq(xp)
            xp = xp[sptr_p]
            colcount_p = n_elements(xp) 
            ymax_col = fltarr(colcount_p)

            for z =0,colcount_p-1 do begin 

             ymax_col[z] = max(image(xp[z],*))

            endfor           

            ymax_col = ymax_col(where(finite(ymax_col))) 
            pthres = max(ymax_col)/ft

            print,'pthres =', pthres
            wait,3

          endif 
 
        endif 

       endif


; Remove old flagging

        ;ptr = where(mask eq pulldownflag, ocount)
	;if (ocount gt 0) then mask[ptr] = mask[ptr] - pulldownflag
	;sxaddhist, name + ' Pulldown bit reset', hm


; Now reflag the pulldown


; First identify all of the triggering pixels
;	ptr = where(image ge pthres, count)

;	sxaddpar, hm, 'NPULL', colcount_p, ' Number of pulldown triggers'
; If no pixels above pulldown threshold, inform user
        ; COMMENTED BY ROBERTA (03/19/2010)
	if (ocount_n eq 0) or (ocount_p eq 0) then begin  
		mstr = 'No pulldown ' 
                mstr = mstr + ' for ' + strn(ifile) ; EDITED BY ROBERTA (03/22/2010)
                print, mstr
                ;
                pthres = 0 & ft = 0 & backg = 0

	endif

 ;IF (COUNT NE 0) THEN BEGIN 	; EDITED BY ROBERTA (03/19/2010)
 if (ocount_p ne 0) then begin

;   mask the columns affected by pulldown
	mask[xp, *] = mask[xp, *] or pulldownflag
        if (write eq 1) then writefits, 'mask_r_ch1_dce1.fits', mask, hm
; and add notes to the mask header
	for i = 0, colcount_p-1 do begin
;		cstr = string(FORMAT='(I03)', i)
;		sxaddpar, hm, 'PXPOS' + cstr, xp[i], ' X position of trigger pixel'
;		sxaddpar, hm, 'PYPOS' + cstr, yp[i], ' Y position of trigger pixel'			
;		sxaddpar, hm, 'PFLX' + cstr, image[trigger_pos[i]], $
;		                  ' Brightness of trigger'			
	endfor		
	
; Next, we estimate what the image without artifacts should look like
; First generate boolean mask for pixels to interpolate
	bptr = where(mask and bad_pix_bits, bitcount)
	bmask = mask * 0B
	if (bitcount gt 0) then bmask[bptr] = 1B
; Output the boolean mask if desired -- for debugging purposes		
	;if (write eq 1) then writefits, 'bmask_r_ch1_dce1.fits', bmask, hm

; Weights for Gaussian smoothing kernel
	;irac_p_make_gweights, 5, 5, sgauss[ch], w5  ; EDITED BY ROBERTA (03/22/2010)
	;irac_p_make_gweights, 11, 11, sgauss[ch], w11 ; EDITED BY ROBERTA (03/22/2010)
        irac_p_make_gweights, 5, 5, sgauss[ch-1], w5
        irac_p_make_gweights, 11, 11, sgauss[ch-1], w11

; Once the flagging has been performed, fill in the masked pixels using
; some estimation scheme, usually a local interpolation but other schemes
; are possible too

; This is the default method using a bootstrapped Gaussian interpolate 
; and a last pass with a larger kernal
	irac_p_estimate_sky, image, bmask, w5, w11, estimated, 	/BOOTSTRAP, $
	                     FALLBACK_INTERP='NAN
        if (write eq 1) then writefits, 'estimated_r32871168_ch1_dce1.fits', estimated, h

; compute estimnate of background in interpolated image 

        backg = median(estimated)

; To use a median filtering based on mmm
;	sky_background, image, bmask, estimated

; To use mmm filtering using an 11x11 window
;	sky_background, image, bmask, estimated, WINDOW=11

; To use a one pass 11x11 Gaussian kernel
;	irac_p_estimate_sky, estimated, bmask, w5, w11, estimated, /SKIP5

; To use a polynomial interpolation of the image
;	irac_p_estimate_sky, image, bmask, w5, w11, estimated, POLY_ORDER=5, $
;	                     FALLBACK_INTERP='POLY', /FALLBACK_ONLY;, SIG=sigma

; You could tinker with the nearest neighbor code as well if you wanted to get 
; very local
;	use_neighbors_for_mux, image, tmask, estimated 

; Output estimated truth image for test purposes
	if (write eq 1) then writefits, 'etruth.fits', estimated, h
	
; Create difference image should just contain the artifacts, if the estimation 
; is good.
	diff = image - estimated

; Now fit the column pulldown, using the flagging from above, instead of
; re-flagging from the mask
		
; 
	for i = 0,colcount_p-1 do begin
                print,'xp(i) =', xp(i)
; Make a vector of the affected column
		vec = diff[xp[i], *]
                vec_im = image[xp[i], *]
; and sigma and mask vectors
		if (keyword_set(USESIG_WEIGHT)) then svec = sigma[xp[i], *] $
		else svec = vec * 0.0D + 1.0D
		mvec = mask[xp[i], *]
; and a row index
		index = indgen(ny)
		
; Replace bad pixels (all but those column pulldown effected) with NaNs
		vptr = where(mvec and bad_notpull_bits, vcount)
		if (vcount gt 0) then vec[vptr] = !VALUES.D_NAN
; Track good pixels
		gptr = where(vec eq vec, gcount)

; If there are no good pixels left in the column, then need to go to the next 
; one -- cheap trick to use a continue statement	
		if (gcount eq 0) then continue

; Basic statistics of column		
		m0 = median(vec)
		mom = moment(vec, /NAN)
		osig = sqrt(mom[1])
		
; Robust estimation of sigma
		sig = robust_sigma(vec)	
; Uncertainty in median is estimated as uncertainty in mean
		sm0 = sig / sqrt(gcount)
		
; Original correction for pulldown was a constant for the entire column
;		corr_vec = vec * 0. + m0
		
; Find maximum pixel in column, this is the dominant pulldown trigger
		cptr = where(vec eq max(vec, /NAN))
; Locate row of maximum pixel in column		
		ymax = index[cptr[0]]
                print,'ymax =', ymax
		
; Count top 5 pixels as a metric -- not used but could be if we try an apriori 
; solution, heritage code at this point
		sptr = bsort(vec[gptr], /REVERSE)
		top5 = total(vec[gptr[sptr[0:4]]])		
		
; Make separate corrections above and below trigger
; Top of column
		topptr = where(index gt ymax, tcount)
                print,'tcount =', tcount
; Bottom of column
		botptr = where(index lt ymax, bcount)
		
; Initialize correction vector
		corr_vec = vec * 0.0	
; If enough pixel above source exist, then fit a 2nd order polynomial to them, 
; require polynomial order + 2 good data points

; Work on top first
		;tgptr = where(vec[topptr] eq vec[topptr], tgcount)  ; EDITED BY ROBERTA (03/20/2010)
               
                IF (TCOUNT NE 0) THEN BEGIN ; EDITED BY ROBERTA (03/20/2010) 

                
                tgptr = where(vec[topptr] eq vec[topptr], tgcount)
                if (tgcount gt norder+2) then begin
; In addition should probably clip outliers before fitting: apply 3-sigma clipping ; EDITED BY ROBERTA (03/26/2010)

                  ;  First, do basic statistics       ; EDITED BY ROBERTA (03/26/2010) 
                  m0_top = median(vec[topptr[tgptr]]) ; EDITED BY ROBERTA (03/26/2010)
                  mom_top = moment(vec[topptr[tgptr]], /NAN) ; EDITED BY ROBERTA (03/26/2010)
                  osig_top = sqrt(mom[1])                     ; EDITED BY ROBERTA (03/26/2010)
   
                  ; Then, do robust estimation of sigma       ; EDITED BY ROBERTA (03/26/2010)
                  sig_top = robust_sigma(vec[topptr[tgptr]])  ; EDITED BY ROBERTA (03/26/2010)
                  ; Uncertainty in median is estimated as uncertainty in mean  ; EDITED BY ROBERTA (03/26/2010)
                  sm0_top = sig_top / sqrt(tgcount)                            ; EDITED BY ROBERTA (03/26/2010)
                             
                  ; now identify and keep only values within a given sigma        ; EDITED BY ROBERTA (03/26/2010)

                  for w =  0, 30 do begin 

                    if (w eq 0) then top_ct_pos = 0 & im_top_ct_neg = 0
                    if (w gt 0) then begin 

                      top_corr_pos =  where((corr_vec[topptr[tgptr]]) > 0, top_ct_pos) 
                      im_top_corr = where(vec_im[topptr[tgptr]] - corr_vec[topptr[tgptr]] < 0, im_top_ct_neg)

                    endif

                    if (top_ct_pos eq 0) and (im_top_ct_neg eq 0) then begin 

                        if (w lt 28) then begin

                         sig_thr = 31. - w

                        endif else begin 

                          sig_thr = 3

                        endelse 

                        tgptr_3s = where(vec[topptr[tgptr]] gt m0_top-sig_thr*sm0_top and vec[topptr[tgptr]] lt  m0_top+sig_thr*sm0_top, tp3scount)   ; EDITED BY ROBERTA (03/26/2010)
 

; Make initial guess at coefficients
                        mt = median(vec[topptr[tgptr[tgptr_3s]]]) ; EDITED BY ROBERTA (03/26/2010)
			pt0 = [mt, dblarr(norder)]

; place fitting data in structure

			ftt = {FLUX:vec[topptr[tgptr[tgptr_3s]]], X:index[topptr[tgptr[tgptr_3s]]]-ymax, $
			       ERR:svec[topptr[tgptr[tgptr_3s]]]}
; Call mpfit to fit polynomial
                        pt = mpfit('exp_func_for_mpfit', pt0, FUNCTARGS=ftt, /quiet)

; Assuming it returns, apply the fit
                        corr_vec[topptr[tgptr]] = pt(0) + pt(1) * exp(-(index[topptr[tgptr]]-ymax)^pt(2))		
                        ;print,'pt(0), pt(1), pt(2) =', pt(0), pt(1), pt(2)  

                    endif else begin 

                     if (w ge 2) then begin

                      w = w - 2

                     endif else  w = w  -1

                      if (w le 28) then begin 

                       sig_thr = 31. - w 

                      endif else begin 

                        sig_thr = 3.
 
                      endelse
                      tgptr_3s = where(vec[topptr[tgptr]] gt m0_top-sig_thr*sm0_top and vec[topptr[tgptr]] lt  m0_top+sig_thr*sm0_top, tp3scount)   ; EDITED BY ROBERTA (03/26/2010)
 
; Make initial guess at coefficients
                        mt = median(vec[topptr[tgptr[tgptr_3s]]]) ; EDITED BY ROBERTA (03/26/2010)
                        pt0 = [mt, dblarr(norder)]

; place fitting data in structure

                        ftt = {FLUX:vec[topptr[tgptr[tgptr_3s]]], X:index[topptr[tgptr[tgptr_3s]]]-ymax, $
                               ERR:svec[topptr[tgptr[tgptr_3s]]]}
; Call mpfit to fit polynomial
                        pt = mpfit('exp_func_for_mpfit', pt0, FUNCTARGS=ftt, /quiet)

; Assuming it returns, apply the fit
                        corr_vec[topptr[tgptr]] = pt(0) + pt(1) * exp(-(index[topptr[tgptr]]-ymax)^pt(2))
                        ;print,'pt(0), pt(1), pt(2) =', pt(0), pt(1), pt(2)

                        w = 30

                     endelse 

;set_plot,'x'
;window,0
;loadct,12
;plot, index[topptr[tgptr[tgptr_3s]]]-ymax, vec[topptr[tgptr[tgptr_3s]]]
;oplot, index[topptr[tgptr[tgptr_3s]]]-ymax, corr_vec[topptr[tgptr[tgptr_3s]]], linestyle=3, color=200
;legend,['sigma = ' + strtrim(sig_thr,2)],/right,/top, box=0, charsize = 1.3
;wait,2


                   endfor 


                   loadct,12
                   set_plot,'ps'
                   device,filename='top_ch1_dce' + (strtrim(float(sxpar(h, 'EXPID')),2)) + '_col' + (strtrim(xp(i),2)) + '.ps',/color
                   plot, index[topptr[tgptr]]-ymax, vec[topptr[tgptr]], yrange=[-1,1],ysty=1,$
                   xtitle='!3 offset from triggering source [pixel]', ytitle='!3 Flux [MJy/sr]', charsize = 1.3, $
                   title='!3 Warm - column '  + strtrim(xp(i),2)
                   oplot, index[topptr[tgptr]]-ymax, corr_vec[topptr[tgptr]], linestyle=3, color =200
                   legend,['top row - column number ' + strtrim(xp(i),2)],/right,/top, box=0, charsize = 2
                   device,/close


; If not enough good data then default to median of top portion for correction
		;endif else if (tgcount gt 0) then corr_vec[tptr] = mt ; EDITED BY ROBERTA (03/20/2010)
                 endif else begin 
   
                   mt = median(vec[topptr])   ; EDITED BY ROBERTA/SEAN (03/22/2010)
                   if (tgcount gt 0) then begin 

                    print,'I am here - top', xp[i], mt
                    corr_vec[tgptr] = mt
                    pt = fltarr(3)
                    pt(0) = mt
                    pt(1) = 0
                    pt(2) = 0

                   endif
 
                 endelse

		
; Now same procedure on bottom of column	
		;bgptr = where(vec[botptr] eq vec[botptr], bgcount)   ; EDITED BY ROBERTA (03/20/2010)
 
                IF (BCOUNT NE 0) THEN BEGIN ; EDITED BY ROBERTA (03/20/2010)

                bgptr = where(vec[botptr] eq vec[botptr], bgcount)
		;IF (BGCOUNT GT 0) THEN mb = median(vec[botptr[bgptr]]) ;  EDITED BY ROBERTA (03/20/2010)
		if (bgcount gt norder+2) then begin
; In addition should probably clip outliers before fitting: apply 3-sigma clipping ; EDITED BY ROBERTA (03/26/2010)

                  ;  First, do basic statistics            ; EDITED BY ROBERTA (03/26/2010)
                  m0_bot = median(vec[botptr[bgptr]])      ; EDITED BY ROBERTA (03/26/2010)
                  mom_bot = moment(vec[botptr[bgptr]], /NAN)   ; EDITED BY ROBERTA (03/26/2010)
                  osig_bot = sqrt(mom_bot[1])                  ; EDITED BY ROBERTA (03/26/2010)
   
                  ; Then, do robust estimation of sigma        ; EDITED BY ROBERTA (03/26/2010)
                  sig_bot = robust_sigma(vec[botptr[bgptr]])    ; EDITED BY ROBERTA (03/26/2010)
                  ; Uncertainty in median is estimated as uncertainty in mean   ; EDITED BY ROBERTA (03/26/2010)
                  sm0_bot = sig_bot / sqrt(bgcount)                             ; EDITED BY ROBERTA (03/26/2010)

                  ; now keep only values > or < a given sigma 

                  for w =  0, 30 do begin

                    if (w eq 0) then bot_ct_pos = 0 & im_bot_ct_neg = 0
                    if (w gt 0) then begin

                      bot_corr_pos =  where((corr_vec[botptr[bgptr]]) > 0, bot_ct_pos)
                      im_bot_corr = where(vec_im[botptr[bgptr]] - corr_vec[botptr[bgptr]] < 0, im_bot_ct_neg)

                    endif

                    if (bot_ct_pos eq 0) and (im_bot_ct_neg eq 0) then begin

                       if (w le 26) then begin 

                        sig_thr = 31. - w

                       endif else begin 

                        sig_thr = 5.

                       endelse 
                       bgptr_3s = where(vec[botptr[bgptr]] gt m0_bot-sig_thr*sm0_bot and vec[botptr[bgptr]] lt  m0_bot + sig_thr*sm0_bot, bt3scount)   ; EDITED BY ROBERTA (03/26/2010)

; Make initial guess at coefficients 
                        mb = median(vec[botptr[bgptr[bgptr_3s]]]) ;  EDITED BY ROBERTA (03/20/2010)
			pb0 = [mb, dblarr(norder)]



; place fitting data in structure
			fb = {FLUX:vec[botptr[bgptr[bgptr_3s]]], X:ymax-index[botptr[bgptr[bgptr_3s]]], $
			       ERR:svec[botptr[bgptr[bgptr_3s]]]}
; Call mpfit to fit polynomial

                        pb = mpfit('exp_func_for_mpfit', pb0, FUNCTARGS=fb,/quiet)  

; Assuming it returns, apply the fit
                        corr_vec[botptr[bgptr]] = pb(0) + pb(1) * exp(-(ymax-index[botptr[bgptr]])^pb(2))

                   endif else begin

                    if (w ge 2) then begin

                      w = w - 2

                    endif else w = w -1

                      if (w le 26) then begin

                       sig_thr = 31. - w

                      endif else begin

                       sig_thr = 5.

                      endelse 
                      bgptr_3s = where(vec[botptr[bgptr]] gt m0_bot-sig_thr*sm0_bot and vec[botptr[bgptr]] lt  m0_bot + sig_thr*sm0_bot, bt3scount)   ; EDITED BY ROBERTA (03/26/2010)

; Make initial guess at coefficients
                        mb = median(vec[botptr[bgptr[bgptr_3s]]]) ;  EDITED BY ROBERTA (03/20/2010)
                        pb0 = [mb, dblarr(norder)]


; place fitting data in structure
                        fb = {FLUX:vec[botptr[bgptr[bgptr_3s]]], X:ymax-index[botptr[bgptr[bgptr_3s]]], $
                               ERR:svec[botptr[bgptr[bgptr_3s]]]}


; Call mpfit to fit polynomial

                        pb = mpfit('exp_func_for_mpfit', pb0, FUNCTARGS=fb,/quiet)   ; EDITED BY ROBERTA (04/12/2010)

; Assuming it returns, apply the fit
                        corr_vec[botptr[bgptr]] = pb(0) + pb(1) * exp(-(ymax-index[botptr[bgptr]])^pb(2))

                        w = 30

                    endelse

;set_plot,'x'
;window,0
;loadct,12
;plot, index[botptr[bgptr[bgptr_3s]]]-ymax, vec[botptr[bgptr[bgptr_3s]]]
;oplot, index[botptr[bgptr[bgptr_3s]]]-ymax, corr_vec[botptr[bgptr[bgptr_3s]]], linestyle=3, color=200
;legend,['sigma = ' + strtrim(sig_thr,2)],/right,/top, box=0, charsize = 1.3
;wait,2

                  endfor
                        
                        
;                 loadct,12
;                 set_plot,'ps'
;                 device,filename='bottom_ch1_dce' + (strtrim(float(sxpar(h, 'EXPID')),2)) + '_col' + (strtrim(xp(i),2)) + '.ps',/color
;                 plot, index[botptr[bgptr]]-ymax, vec[botptr[bgptr]], yrange=[-1,1],ysty=1, $
;                 xtitle='!3 offset from triggering source [pixel]', ytitle='!3 Flux [MJy/sr]', charsize = 1.3, $
;                 title='!3 Warm - column ' + strtrim(xp(i),2)
;                 oplot, index[botptr[bgptr]]-ymax, corr_vec[botptr[bgptr]], linestyle=3, color=200		
;                 legend,['bottom row - column number ' + strtrim(xp(i),2)],/right,/top, box=0, charsize = 2 
;                 device,/close
                         

; If not enough good data then default to median of top portion for correction
		endif else begin 
            
                  mb = median(vec[botptr]) ;  EDITED BY ROBERTA (03/20/2010)
                  print,'bgcount =', bgcount
                  print,'vec[botptr] =',vec[botptr]
                  if (bgcount gt 0) then begin 

                   print,'I am here - bottom', xp[i], mb
                   corr_vec[bgptr] = mb
                   pb = fltarr(3)
                   pb(0) = mb
                   pb(1) = 0
                   pb(2) = 0

                  endif

                endelse

;                  loadct,12
;                  set_plot,'ps'
;                  device,filename='top_bottom_ch1_dce' + (strtrim(float(sxpar(h, 'EXPID')),2)) + '_col' + (strtrim(xp(i),2)) + '.ps',/color
;                  plot, findgen(255)-ymax, vec, xrange=[min(findgen(255)-ymax),170], xsty=1 ,yrange=[-1,1],ysty=1, $
;                  xtitle='!3 offset from triggering source [pixel]', ytitle='!3 Flux [MJy/sr]', charsize = 1.3, charthick=1.5, $
;                  title='!3 Warm - column '  + strtrim(xp(i),2)
;                  oplot, findgen(255)-ymax,  corr_vec, linestyle=3, color = 200
;                   legend,['F!Ibottom!N = ' + strtrim(pb(0),2) + ' + ' + strtrim(pb(1),2) + ' * e!Ex^' + strtrim(pb(2),2) + '!N',$
;                  '', 'F!Itop!N = ' + strtrim(pt(0),2) + ' + ' + strtrim(pt(1),2) + ' * e!Ex^' + strtrim(pt(2),2) + '!N'],/right,/top, box=0, charsize = 1.3
;                  device,/close
                

              ENDIF  ; EDITED BY ROBERTA (03/19/2010)


		
; Now average adjacent rows' solutions to fit trigger pixel	
		if (ymax eq 0) then corr_vec[ymax] = corr_vec[ymax + 1] else $
		if (ymax eq ny-1) then corr_vec[ymax] = corr_vec[ymax - 1] else $ 	
		corr_vec[ymax] = (corr_vec[ymax + 1] + corr_vec[ymax - 1]) / 2.0
		
; Update image with correction


		image[xp[i], *] = image[xp[i], *] - corr_vec
; and update sigma
		sigma[xp[i], *] = sigma[xp[i], *] * pullnoise_factor
		
; and spew to log file
;		if (keyword_set(PLOG)) then printf, plun, xt[i], $
;			                              ymax, cmax, top5, m0, sig, mtop, $
;										  smtop, mbot, smbot, a1, b1, a2, b2


          ENDIF  ; EDITED BY ROBERTA (03/19/2010)

	endfor	

; Update header 
	sxaddpar, h, 'PULLCOR', 'T', ' Pulldown corrected'
	sxaddpar, hs, 'PULLCOR', 'T', ' Pulldown corrected'
; 
	nmstr = ' Noise multiplier for pulldown correction'
	sxaddpar, hs, 'PULLSIG', pullnoise_factor, nmstr

; Output corrected image, uncertainty and updated mask
	writefits, ofile, image, h
	writefits, osfile, sigma, hs
	writefits, omfile, mask, hm
	



 ENDIF ELSE BEGIN ; EDITED BY ROBERTA (03/19/2010) 

        writefits, ofile, image, h
        writefits, osfile, sigma, hs
        writefits, omfile, mask, hm

 ENDELSE

; If open, close log files -- not currently working  ; EDITED BY ROBERTA (03/26/2010)
        if (keyword_set(PLOG)) then free_lun, plun   ; EDITED BY ROBERTA (03/26/2010)

endfor
return  ; EDITED BY ROBERTA (03/26/2010)

end

pro segment_list, list, nsegments, seg_start, seg_stop, DISTANCE=distance, $
                  WRAP_VALUE=wrap_value, START_INDEX=start_index, $
				  STOP_INDEX=stop_index
;+
; NAME:
;     SEGMENT_LIST
;
; PURPOSE:
;      Returns number of and start and stop positions of contiguous segments 
;      in a list of pixels. 
;
; SYNTAX:
;      SEGMENT_LIST, list, nsegments, seg_start, seg_stop, DISTANCE=distance, $
;                    WRAP_VALUE=wrap_value, START_INDEX=start_index, $
;                    STOP_INDEX=stop_index
; INPUT:
;      list: array containing positions
; OUTPUTS:
;      nsegments: number of independent segments in list
;      seg_start: array containing values of segment starts 
;      seg_stop: array containing values of segment end points
;
; OPTIONAL INPUT KEYWORDS:
;      DISTANCE: maximum distance in index units for points to be considered
;                contiguous
;      WRAP_VALUE: wrap around value
;
; OPTIONAL OUTPUT KEYWORDS:
;      START_INDEX: indices of segment start points
;      STOP_INDEX: indices of segment stop points
;
; NOTES:
;      Need to handle array wraparound
; HISTORY:
;      Initial code doodle 31 Mar 2005 SJC
;      Properly handle case of only one segment 11 April 2005 SJC
;-

; Default separation
	if (not keyword_set(DISTANCE)) then distance = 1
	nlist = n_elements(list)
; Find spacing between adjacent elements in list
	sep = list - shift(list, +1)

; Check for segments, not including wraparound
	ptr = where(sep gt distance, count)
	nsegments = count + 1

; Now set the beginning points of the segments
	if (count gt 0) then begin
		start_index = [0, ptr]
		seg_start = [list[0], list[ptr]]

; and end segments
		stop_index = [ptr-1L, nlist-1L]
		seg_stop = [list[ptr-1L], list[nlist-1L]]

	endif else begin
; entire index is one segment
		start_index = [0]
		seg_start = [list[0]]

; and end segments
		stop_index = [n_elements(list)-1L]
		seg_stop = [list[n_elements(list)-1L]]
	endelse
	
; If wrap value is set, see if we do indeed wrap around
	if (keyword_set(WRAP_VALUE) and nsegments gt 1) then begin
		xsep0 = sep[0] + wrap_value
; Is the first entry contiguous with wrapped data
		if (xsep0 le distance) then begin
; then the number of segments is one less than calculated
; the first list entry is not an endpoint, but the endpoint of this 
; segment is really the endpoint of the very last segment
; will need to check output to see if there are any endpoints that are less
; than the start points
			seg_start = seg_start[1L:nsegments-1L]
			start_index = start_index[1L:nsegments-1L]
			seg_stop = [seg_stop[1:nsegments-2L], seg_stop[0]]
			stop_index = [stop_index[1:nsegments-2L], stop_index[0]]
; Decrement number of segments after merging first and last
			nsegments = nsegments - 1
		endif
	endif

return
end

pro irac_p_estimate_sky, image, mask, w5, w11, iest, BOOTSTRAP=bootstrap, $
                         FALLBACK_INTERP=fallback_interp,  MAXITER5=maxiter5, $
                         MAXITER11=maxiter11, SKIP5=skip5, $
                         FALLBACK_ONLY=fallback_only, POLY_ORDER=poly_order, $
                         SIG=sig
;+
; PROCEDURE:
;   IRAC_P_ESTIMATE_SKY
;
; SYNTAX:
;   IRAC_P_ESTIMATE_SKY, image, mask, w5, w11, image_estimated [, /BOOTSTRAP, $
;                        FALLBACK_INTERP={}, MAXITER5=maxiter5,$
;                        MAXITER11=maxiter11, /SKIP5, /FALLBACK_ONLY]
;
; INPUTS:
;    image: image array
;    mask: boolean mask for image 
;    w5: 5x5 pixel smoothing kernel
;    w11: 11x11 pixel smoothing kernel
; 
; OUTPUT:
;    iest: estimated sky image
; 
; KEYWORDS:
;     BOOTSTRAP:
;     FALLBACK_INTERP
;     MAXITER5
;     MAXITER11
;     SKIP5 : If set, do not perform the 5x5 pixel interpolation
;     FALLBACK_ONLY: If set, perform only the fallback interpolation
;
; PURPOSE:
;    To estimate the sky for pixels which have been masked off.  The
;    resulting estimate should be good enough to perform artifact fitting 
;    through the data.
;
; METHOD:
;    For each masked pixel, find the number of unmasked good pixels in a 5x5
;    box around the pixel.  If the number of pixels is greater than a
;    threshold, replace the masked pixel value with the weighted average of 
;    the good unmasked pixels in the box.  The weight image is passed as an
;    input to the procedure and is usually a Gaussian weight.  If the 5x5 box
;    does not contain enough good unmasked pixels, the same procedure is
;    applied using an 11x11 box.  If the 11x11 box does not contain enough
;    good unmasked pixels, then the masked pixel value is determined by the 
;    fallback interpolation method specified by the FALLBACK_INTERP keyword.
;    
;    If no fallback option is specified, then the pixel value is replaced with
;    a NaN.  The current fallback interpolation options are to replace the
;    pixel with the median of the image and to replace the pixel with an
;    interpolated value using the tri_surf procedure.
;
; FUNCTIONS USED:
;    TRI_SURF(), MEDIAN(), TOTAL()
;
; HISTORY: 
;    07 Sep 2004 SJC Initial headscratching
;    22 Nov 2004 SJC 
;    23 Mar 2005 SJC Reduced amount of data needed to apply weighting
;    07 June 2005 SJC Correct bug that used channel 4 weights for channel 3
;                  data (WTR found)
;    12 Aug 2005 SJC Removed from main loop
;    19 Aug 2005 SJC Added FALLBACK_INTERP keyword which controls now to 
;                    interpolate the data if there is not enough info
;    10 Mar 2006 SJC Added BOOTSTRAP keyword, if there are pixels that are
;                    for which a background is not initially estimated, the 
;                    previous pass is smoothed to determine the background 
;                    until all pixels are estimated
;    19 Feb 2008 SJC Revised to properly implement the bootstrap method. Fairly ;                    significant code changes
;    25 Nov 2008 SJC Added options to perform a less local estimation of the sky
;                
;     
;-

	name = 'IRAC_P_ESTIMATE_SKY:'

; Parse parameters
	if (N_params() ne 5) then begin
		print,name+': Syntax -- IRAC_P_ESTIMATE_SKY, image, mask, w5, w11, $'
		print,name+':                                iest, /BOOTSTRAP, $ '
		print,name+':      FALLBACK_INTERP=[NAN, TRI_SURF, MEDIAN, POLY, MEAN]'
		print,name+': default fallback interp is to replace masked pixels with '
		print,name+': image mean'
; GOO Need to add a little more documentation here
		return
	endif

; Default fallback interpolation scheme is not to interpolate at all
	if (not keyword_set(FALLBACK_INTERP)) then fallback_interp = 'MEAN'

; Need control option if bootstrap is not applied
	if (keyword_set(BOOTSTRAP)) then do_bootstrap = 1 else do_bootstrap = 0
	
; Maximum number of iterations for bootstrap method
	if (not keyword_set(MAXITER5)) then maxiter5 = 20
	if (not keyword_set(MAXITER11)) then maxiter11 = 2
	if (not keyword_set(SKIP5)) then skip5 = 0
	if (not keyword_set(FALLBACK_ONLY)) then fallback_only = 0

; Minimum number of valid pixels to use a specified smoothing window
	min_good_5by5 = 5
	min_good_11by11 = 11

; NaN value
	nan = !VALUES.F_NAN

; Image dimensions
	sz = size(image)
	nx = sz[1]
	ny = sz[2]

; Create intrinsic weight file -- 
	iw = lonarr(nx, ny)
	gptr = where(finite(image) and mask eq 0B, count)
	
; If there are good pixels in the image, then let's try to make an estimate of
; what should be in the masked pixels
	if (count gt 0) then begin
; set Boolean weight for input image
		iw[gptr] = 1L
	endif else begin
; No unmasked good pixels exist in image -- estimation cannot work
		print, name+': No good data exists in image -- returning a NaNed image'
		iest = dblarr(nx, ny) + nan
		return
	endelse

; Estimated image is the original image with masked pixels replaced, start off 
; with the original image
	iest = image

; Loop over all masked pixels, for each pixel find the pixels in the
; adjacent area and replace the affected pixel with the weighted sum of 
; the adjacent pixels, continue the loop until we all pixels are estimated or
; a maximum number of iterations is exceeded.  
	done = 0
	iter5 = 0L
	
	while (not done and not skip5 and not fallback_only) do begin
; Find all pixels to estimate pixels
		ptr = where(iw eq 0B, count)
; and their locations on the array
		x = ptr mod nx
		y = ptr / ny
		
; Create a temporary weight array which holds the current bad pixel state for 
; this iteration throught the image 
		tw = iw
; and a temporary estimate image which holds the initial state of the image for 
; this pass
		tiest = iest
		
; Intialize counter for pixels to estimate
		i = 0L
; and perform a pass over all unestimated pixels
		while (i le count-1) do begin

; Start with a 5x5 pixel region, trim if at edges of array
			xa = (x[i] - 2) > 0
			xb = (x[i] + 2) < (nx - 1L)
			ya = (y[i] - 2) > 0
			yb = (y[i] + 2) < (ny - 1L)
; array elements of corresponding portion of smoothing kernel (which is centered
; on pixel position (2,2)
			ua = xa - x[i] + 2
			ub = xb - x[i] + 2
			va = ya - y[i] + 2
			vb = yb - y[i] + 2

; See if there are enough pixels 
			ngood_pixels = total(tw[xa:xb, ya:yb])

			if (ngood_pixels gt min_good_5by5) then begin
; Grab subimage, sub kernel, need to ignore any NaNs in image
				siw = tw[xa:xb, ya:yb]
				sw5 =  w5[ua:ub, va:vb]
				siest = tiest[xa:xb, ya:yb]
; Need to ignore NaNs in the image
				gptr = where(siest eq siest, gcount)
; Now replace the pixel with the kernel convolved average of the surrounding 
; good pixels
				if (gcount gt 0) then begin
					subweight = siw[gptr] * sw5[gptr]
					iest[ptr[i]] = total(siest[gptr] * subweight) $
			         	                / total(subweight)
; and update the weight for the pixel as it is good for the purposes of the next ; iteration
					iw[ptr[i]] = 1B
				endif 
			endif
; go to the next pixel in the list
			i = i + 1L
		endwhile
; Now check to see if we need to perform another pass, increment counter and see ; if any bad pixels are remaining
		iter5 = iter5 + 1
		print, 'iter5 = ', iter5
		ptr = where(iw eq 0B, count)
		if (iter5 ge maxiter5 or count eq 0 or do_bootstrap eq 0) then done = 1
	endwhile

; If there are any bad pixels remaining, try to fill them in using a larger
; smoothing kernel.
; Loop over all masked pixels, for each pixel find the pixels in the
; adjacent area and replace the affected pixel with the weighted sum of 
; the adjacent pixels, continue the loop until we all pixels are estimated or
; a maximum number of iterations is exceeded.  
	if (count gt 0) then done = 0 else done = 1
	iter11 = 0L
		
	while (not done and not fallback_only) do begin
; Find all pixels to estimate pixels
		ptr = where(iw eq 0B, count)
; and their locations on the array
		x = ptr mod nx
		y = ptr / ny
		
; Create a temporary weight array which holds the current bad pixel state for 
; this iteration throught the image 
		tw = iw
; and a temporary estimate image which holds the initial state of the image for 
; this pass
		tiest = iest
		
; Intialize counter for pixels to estimate
		i = 0L
; and perform a pass over all unestimated pixels
		while (i le count-1) do begin

; Start with a 11x11 pixel region, trim if at edges of array
			xa = (x[i] - 5) > 0
			xb = (x[i] + 5) < (nx - 1L)
			ya = (y[i] - 5) > 0
			yb = (y[i] + 5) < (ny - 1L)
; array elements of corresponding portion of smoothing kernel (which is centered
; on pixel position (5,5)
			ua = xa - x[i] + 5
			ub = xb - x[i] + 5
			va = ya - y[i] + 5
			vb = yb - y[i] + 5

; See if there are enough pixels 
			ngood_pixels = total(tw[xa:xb, ya:yb])

			if (ngood_pixels gt min_good_11by11) then begin
; Grab subimage, sub kernel, need to ignore any NaNs in image
				siw = tw[xa:xb, ya:yb]
				sw11 =  w11[ua:ub, va:vb]
				siest = tiest[xa:xb, ya:yb]
; Need to ignore NaNs in the image
				gptr = where(siest eq siest, gcount)
; Now replace the pixel with the kernel convolved average of the surrounding 
; good pixels
				if (gcount gt 0) then begin
					subweight = siw[gptr] * sw11[gptr]
					iest[ptr[i]] = total(siest[gptr] * subweight) $
			         	                / total(subweight)
; and update the weight for the pixel as it is good for the purposes of the next ; iteration
					iw[ptr[i]] = 1B
				endif 
			endif
; go to the next pixel in the list
			i = i + 1L
		endwhile
; Now check to see if we need to perform another pass, increment counter and see ; if any bad pixels are remaining
		iter11 = iter11 + 1
		ptr = where(iw eq 0B, count)
		print, 'iter11 ', iter11
		if (iter11 ge maxiter11 or count eq 0 or do_bootstrap eq 0) then done=1

	endwhile
	
; Moved the fallback estimation to the end of the scheme.  Assuming the 
; estimated image isn't too bad, let's use that in the fallback scheme
		
; If the weighted interpolation window scheme fails due to a lack of
; good data, a fallback needs to be used.  Currently, there are 3 options,
; interpolate the entire image over the bad pixels using tri_surf, replace the
; bad pixels with the array wide median or replace the bad pixels with NaNs.

; Only continue if bad pixels still exist
	if (count gt 0) then begin
; Start with a vector of all the good pixels
		gptr = where(iw eq 1B, gcount)
		if (fallback_interp eq 'TRI_SURF') then begin
; As a last resort create an interpolated surface and uses pixels from that 
; to create estimated image, tri_surf is memory intensive, may want to remove
; the last resort method needs to be improved
			gx = ptr mod nx
			gy = ptr / ny
			xest = tri_surf(iest[gptr], gx, gy, BOUNDS=[0, 0, nx-1, ny-1], $
		   	             GS=[1,1], /LINEAR)
		endif else if (fallback_interp eq 'MEDIAN') then begin
			xest = dblarr(nx, ny) + median(iest[gptr])
		endif else if (fallback_interp eq 'MEAN') then begin
			xest = dblarr(nx, ny) + total(iest[gptr]) / gcount
		endif else if (fallback_interp eq 'NAN') then begin
			xest = dblarr(nx, ny) + nan
		endif else if (fallback_interp eq 'POLY') then begin
; Set poly surface order		
			if (not keyword_set(POLY_ORDER)) then poly_order = 2
; Set initial guess for surface coefficients
print, 'Poly-fit'
			fname = 'poly2d_func_for_mpfit'
			ncoeff = (poly_order + 1) * (poly_order + 2) / 2
			coefs0 = dblarr(ncoeff)
			coefs0[0] = median(iest[gptr])
; Create array indices
			xx = findgen(nx) # replicate(1.0, ny)
			yy = replicate(1.0, nx) # findgen(ny)
; If uncertainty image exists use it, otherwise give all pixels uniform weight
			if (not keyword_set(SIG)) then sig = iest * 0. + 1.
; Create structure to pass function arguments to mpfit
			fa = {FLUX: iest[gptr], ERR: sig[gptr], X: xx[gptr], Y: yy[gptr]}
			stop
; Fit polynomial surface to good pixels
			coefs = mpfit(fname, coefs0, FUNCTARGS=fa)
; Return surface
			xest = poly2d_func(coefs, X=xx, Y=yy)
; Set default to mean for the moment
		endif else begin
			xest = dblarr(nx, ny) + total(iest[gptr]) / gcount
		endelse
	endif
	
; Now replace the unestimated pixels with this more global estimator
	ptr = where(iw eq 0B, count)
	if (count gt 0) then iest[ptr] = xest[ptr]

return
end

function poly1d_func_for_mpfit, p, FLUX=flux, ERR=err, X=x
; The order of coefficients, p, is the same for this function as
; the IDL function, POLY, p[0] + p[1] * x+ p[2] * x * x + ....
	if (not keyword_set(FLUX) or not keyword_set(X) or $
	    not keyword_set(ERR)) then $
	     message, 'INPUT keywords, FLUX, ERR, X and Y need to be passed'

	nf = n_elements(flux)
	nx = n_elements(x)
	nerr = n_elements(err)
	if (nf ne nx or nf ne nerr) then $
		message, 'Input f, x, err need to have # of elements'

	n = n_elements(p)
	order = n - 1

	model = 0. * flux
	index = 0
	for i = 0, order do begin
		model = model + p[index] * x^i
		index = index + 1
	endfor
	
	model = (flux - model) / err
return, model
end

; EDITED BY ROBERTA 04/12/2010

function exp_func_for_mpfit, p, FLUX=flux, ERR=err, X=x

        if (not keyword_set(FLUX) or not keyword_set(X) or $
            not keyword_set(ERR)) then $
             message, 'INPUT keywords, FLUX, ERR, X and Y need to be passed'

        nf = n_elements(flux)
        nx = n_elements(x)
        nerr = n_elements(err)
        if (nf ne nx or nf ne nerr) then $
                message, 'Input f, x, err need to have # of elements'

        n = n_elements(p)
        model = 0. * flux
        model = model + p[0] + p[1]* exp(-x^p[2])

        model = (flux - model) / err
return, model
end


function poly2d_func_for_mpfit, p, FLUX=flux, ERR=err, X=x, Y=Y
	if (not keyword_set(FLUX) or not keyword_set(X) or $
	     not keyword_set(Y) or not keyword_set(ERR)) then $
	     message, 'INPUT keywords, FLUX, ERR, X and Y need to be passed'

	nf = n_elements(flux)
	nx = n_elements(x)
	ny = n_elements(y)
	nerr = n_elements(err)
	if (nf ne nx or nf ne ny or nf ne nerr) then $
		message, 'Input f, x, y, err need to have # of elements'

	n = n_elements(p)
	order = (-1. + sqrt(1. + 8. * n)) / 2. - 1
	if (order ne fix(order)) then $
		message, $
		'Number of parameters is not appropriate for 2d polynomial expansion'
	
	model = 0. * flux
	index = 0
	for i = 0, order do begin
		for j = 0, i do begin
			model = model + p[index] * x^(i-j) * y^j
;			print, 'index = ', i, ' xpow = ', (i-j), ' ypow =', j, ' order =', i
			index = index + 1
		endfor
	endfor
	
	model = (model - flux) / err
return, model
end

function poly2d_func, p, X=x, Y=Y
	if (not keyword_set(X) or not keyword_set(Y)) then $
	     message, 'INPUT keywords, FLUX, X and Y need to be passed'
		
	nx = n_elements(x)
	ny = n_elements(y)
	
	if (nx ne ny) then $
		message, 'Input x, y need to have # of elements'

	n = n_elements(p)
	order = (-1. + sqrt(1. + 8. * n)) / 2. - 1
	if (order ne fix(order)) then $
		message, $
		'Number of parameters is not appropriate for 2d polynomial expansion'

; Initialize function	
	model = 0. * x
	index = 0
; Build polynomial	
	for i = 0, order do begin
		for j = 0, i do begin
			model = model + p[index] * x^(i-j) * y^j
			index = index + 1
		endfor
	endfor
	
return, model
end

pro irac_p_make_gweights, nx, ny, sigma, weights
;+
; PROCEDURE:
;    IRAC_P_MAKE_GWEIGHTS
;
; SYNTAX:
;    IRAC_P_MAKE_GWEIGHTS, NX, NY, SIGMA, OUTFILE
;
; PURPOSE: 
;    Generate weight files for use in Gaussian interpolation of images.
;    Primarily used in IRAC_P_MAKE_ESTIMATED_IMAGES
;
; HISTORY:
;    Written pre 11 Aug 2005 SJC
;    Adapted to just return matrix of weights from 
;          IRAC_P_MAKE_GAUSSIAN_WEIGHTS07 Oct 09 SJC
;  
;-

; Parse parameters
	if (N_params() ne 4) then begin
		print, 'Syntax -- IRAC_P_MAKE_GWEIGHTS, nx, ny, sigma, outfile'
		return
	endif

; Make x and y dimension odd
	if (nx/2 eq float(nx)/2.0) then nx = nx + 1
	if (ny/2 eq float(ny)/2.0) then ny = ny + 1

; Determine center pixel
	cx = (nx - 1) / 2
	cy = (ny - 1) / 2

; Create weights
	x = (indgen(nx) - cx) # replicate(1, ny)
	y = replicate(1, nx) # (indgen(ny) - cy)
	d = x * x + y * y

; Allocate weight array
	weights = exp(float(-d) / (2.0 * sigma))


return
end


function muxmodel_for_mpfit, params, TRIGGERS=triggers, INDEX=index, $
         DATA=data, UNC=unc
;+
; FUNCTION:
;    MUXMODEL_FOR_MPFIT
;
; SYNTAX:
;    vals = muxmodel_for_mpfit(params, TRIGGERS=triggers, INDEX=index, 
;                               DATA=data, UNC=unc)
;
; PURPOSE: 
;    Model of muxbleed used in conjunction with mpfit procedure to fit 
;    instances of muxbleed.
;
; NOTES:
;    Functional form developed by I. Song and implemented in IRAC BCD pipeline.
;    The scaling of the muxbleed is slightly stochastic with triggering pixel 
;    brightness.  In some cases, the muxbleed is oversubtracted and others 
;    undersubtracted.  For oversubtracted cases, it appears that most of the 
;    oversubtraction occurs in the first few pixels, 
;
; HISTORY:
;    Written 08 Oct 2009 SJC
;-

	ntriggers = n_elements(params)
	if (n_elements(triggers) ne ntriggers) then $
	    message, 'number of triggers does not match number of fit coefficients'
; Number of samples to flag (3 rows worth)
	nsamples_to_flag = 192L
; Coefficients for the fit
	coef = [3.1880D, -2.4973D, 1.2010D, -0.2444D]
; Truncation length
	trunc_point = 4

; Output model has same dimensionality as data
	model = data * 0.0D
; Apply no weight to pixels adjacent to triggers
	weight = model + 1.D

; Build model by looping over triggers
	for i = 0, ntriggers-1 do begin
		tmin = triggers[i] + 4L
		tmax = tmin + (nsamples_to_flag - 1L) * 4L
		ptr = where(index ge tmin and index le tmax, count)
		if (count eq 0) then message, 'No flagged pixels for this trigger'
		x = alog10(dindgen(192) + 1)
		val = coef[0] + coef[1] * x + coef[2] * x * x + coef[3] * x * x * x
; For cases where the muxbleed is oversubtracted, then we want to truncate the 
; fix
		if (params[i] lt 0.) then val[trunc_point: *] = 0.0D
		model[ptr] = model[ptr] + params[i] * 10.D^val
; Apply zero weight to position right next to the trigger		
		weight[ptr[0]] = 0.0D
; Also apply zero weight to triggering pixels in between
		xptr = where(index eq triggers[i], xcount)
		if (xcount eq 1) then weight[xptr[0]] = 0.0D
	endfor

; Return deviation for use by MPFIT
return, weight * (data - model) / unc
end

function muxmodel, params, TRIGGERS=triggers, INDEX=index
;+
; FUNCTION:
;    MUXMODELT
;
; SYNTAX:
;    vals = muxmodel(params, TRIGGERS=triggers, INDEX=index)
;
; PURPOSE: 
;    Model of muxbleed to remove from overlapping instances of data.
;
; HISTORY:
;    Writte 08 Oct 2009 SJC
;-

	ntriggers = n_elements(params)
	if (n_elements(triggers) ne ntriggers) then $
	    message, 'number of triggers does not match number of fit coefficients'
; Number of samples to flag (3 rows worth)
	nsamples_to_flag = 192L
; Coefficients for the fit
	coef = [3.1880D, -2.4973D, 1.2010D, -0.2444D]
; Truncation length
	trunc_point = 4

; Output model has same dimensionality as data
	model = index * 0.0D

; Build model by looping over triggers
	for i = 0, ntriggers-1 do begin
		tmin = triggers[i] + 4L
		tmax = tmin + (nsamples_to_flag - 1L) * 4L
		ptr = where(index ge tmin and index le tmax, count)
		if (count eq 0) then message, 'No flagged pixels for this trigger'
		x = alog10(dindgen(192) + 1)
		val = coef[0] + coef[1] * x + coef[2] * x * x + coef[3] * x * x * x
; For cases where the muxbleed is oversubtracted, then we want to truncate the 
; fix
		if (params[i] lt 0.) then val[trunc_point: *] = 0.0D		
		model[ptr] = model[ptr] + params[i] * 10.D^val
	endfor

; Return fit
return, model
end

pro sky_background, im, mask, back, WINDOW=window
;+
; PURPOSE: 
;     Calculate an estimate of sky background using sky.pro and appropriate
;     filter size for each pixel in an image
;-

	if (not keyword_set(WINDOW)) then window = 16

	nim = im
	ptr = where(mask ne 0, count)
	if (count gt 0) then nim[ptr] = !VALUES.D_NAN
	back = im * 0.0D
	
	sz = size(im)
	if (sz[0] ne 2) then message, 'Only works on 2d images'
	nx = sz[1]
	ny = sz[2]
	
	for i = 0, nx-1 do begin
		x0 = (i - window/2) > 0
		x1 = (i + window/2) < (nx-1)
		for j = 0, ny-1 do begin
			y0 = (j - window/2) > 0
			y1 = (j + window/2) < (ny-1)
			sky, nim[x0:x1, y0:y1], smode, ssig, /SILENT, /NAN
			back[i, j] = smode
		endfor
	endfor

return
end

pro use_neighbors_for_mux, image, mask, back, MUXBIT=muxbit
	if (not keyword_set(MUXBIT)) then muxbit = 5
	ptr = where(mask and 2L^muxbit, count)
	back = image
	if (count gt 0) then begin
		sz = size(image)
		nx = sz[1]
		ny = sz[2]
		xx = lindgen(nx) # replicate(1, ny)
		yy = replicate(1, nx) # lindgen(ny)
		
		for i = 0, count-1 do begin
			x = ptr[i] mod nx
			y = ptr[i] / ny
			xlo = (x - 1) > 0
			xhi = (x + 1) < (nx - 1)
			ylo = (y - 1) > 0
			yhi = (y + 1) < (ny - 1)
			sub = image[xlo:xhi, ylo:yhi]
			gptr = where(sub eq sub, gcount)
			if (gcount gt 0) then image[x, y] = total(sub[gptr]) / (count-1.) $
			                                                    - image[x, y]
		endfor		
	endif
return
end

pro remove_leading_path, istr, ostr, LPATH=lpath
;+
; NAME:
;    REMOVE_LEADING_PATH
;
; PURPOSE:
;    Removes the leading directories from an input string
; 
; INPUT:
;    istr: string containing input pathname
;
; OUTPUT:
;    ostr: string containing output filename
; 
; OUTPUT Keyword
;    lpath: if set, returns the leading path of the input string
;
; METHOD:
;    Find the last '/' in the string and remove all characters from the
;    first to that position.  For lpath, returns from the first character to
;    the last '/'.  If there are no '/', then lpath is a null string.
;
; HISTORY
;    Initial coding SJC 29 Dec 2006
;    Module checked SJC 29 Dec 2006
;    Added leading path option 23 April 2008
;-

; Check input parameters
	if (N_params() ne 2) then begin
		print, 'Syntax -- remove_leading_path, in_string, out_string, $'
		print, '              LPATH=leading_path'
		return
	endif

; Set slash string
    slash_str = '/'

; Copy the input string to the output string
    ostr = istr

; Find the position of the last slash in the input string
    pos = strpos(ostr, slash_str, /REVERSE_SEARCH)

; If the lash slash is not the last character then add a slash to the end of 
; the string
    if (pos ne -1) then begin
        lpath = strmid(ostr, 0, pos+1)
    	ostr = strmid(ostr, pos+1, strlen(ostr)-pos-1) 
    endif else lpath = ''
    
return
end
