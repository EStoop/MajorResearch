pro irac_warm_correct_pulldown_v2, ifile, sfile, mfile, ofile, osfile, $
                             omfile, pthres, ft, backg, PTHRESH=pthresh, $
                             PULLNOISE_FACTOR=pullnoise_factor, $
                             VERBOSE=verbose, USESIG_WEIGHT=usesig_weight
;+
; NAME:
;   IRAC_WARM_CORRECT_PULLDOWN_V2
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
;    PTHRES : if set, the code computes the threshold in MJy/sr * second 
;             for pulldown to occur.
;    PULLNOISE_FACTOR : Floating point, if set this is the value to inflate the
;                       uncertainty image by after pulldown correction
;    USESIG_WEIGHT : Boolean, if set, then use the sigma image for weighting 
;                    the data in polynomial fitting the pulldown
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
;   EXP_FUNC_FOR_MPFIT, POLY
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
;    2 August 2010 Roberta Paladini - new flagging/ fitting technique (exponential function)
;-


; Define constants
	name = 'IRAC_WARM_CORRECT_PULLDOWN' 
	syntax =['irac_warm_correct_pulldown, image, uncertainty, mask,$', $
	         '                                         oimage, ounc, omask, pull_th, int_time, backgr']
	
; SCLK time of warm mission start
	warm_start_sclk = 927192851.D
	
; Threshold for muxbleed to start in MJy * s / sr
        pulldown_thresh = [400, 100]   
; Bit set in imask for muxbleed
	pulldownflag = 2^7

; Width for Gaussian smoothing kernel (sigma = FWHM / 1.665) in pixels
	sgauss = [1.44, 1.43] / 1.665 / 1.22
; Bits set for bad pixels
	bad_pix_bits = long(total(2L^[3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]))
; Set bits for data not to have solution fit to, assume that straylight, 
; muxbleed, banding are alright for this case (so 3, 5, 6 are okay)
	bad_notpull_bits = long(total(2L^[8, 9, 10, 11, 12, 13, 14, 15]))

 	
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

; Set minimum number of points to fit artifact and apply correction
  npoints = 4
	
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


; Read in image file
	image = readfits(ifile, h, SILENT=silent)
	if (n_elements(image) eq 1) then message, 'Not a valid BCD file'
; Should test to make sure it is a BCD, otherwise this procedure is not useful
	naxis = sxpar(h, 'NAXIS')
	if (naxis eq 3) then message, 'Code does not work for subarray BCDs'

; Obtain information about image from header
	ch = sxpar(h, 'CHNLNUM') 
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

	
; Determine threshold (MJy/sr) for pulldown trigger: flagging section of the code/ redefinition of bimsk

        if (keyword_set(PTHRESH)) then begin

        ; first compute robust median of image columns

        super_im_med_col_vec = fltarr(nx)
        for t = 0,nx-1 do begin 

          im_col = image(t,*)

          ;;;;;;;

           medclip, im_col, super_im_med_col, super_im_sig_col, clipsig=3.5, maxiter=5
           super_im_med_col_vec(t) = super_im_med_col

        endfor

        ;;;;;;;; now compute median of all columns: super median

        medclip, super_im_med_col_vec, super_im_med_all, super_im_sig_all, clipsig=3.5, maxiter=5

        sigma_to_med = abs(super_im_sig_all/super_im_med_all)
        
        case 1 of

           (sigma_to_med lt 0.1): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (5*super_im_sig_all)), ocount_p)
           (sigma_to_med ge 0.1) and (sigma_to_med lt 0.3): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (3.*super_im_sig_all)), ocount_p)
           (sigma_to_med ge 0.3) and (sigma_to_med lt 0.5): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (2.*super_im_sig_all)), ocount_p)
           (sigma_to_med ge 0.5) and (sigma_to_med lt 1) and (super_im_med_all ge 0.1): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (2*super_im_sig_all)), ocount_p)
           (sigma_to_med ge 0.5) and (sigma_to_med lt 1) and (super_im_med_all lt 0.1): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (0.5*super_im_sig_all)), ocount_p)
           (sigma_to_med gt 1.) and (sigma_to_med lt 5) and (super_im_med_all ge 0.1): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (0.5*super_im_sig_all)), ocount_p)
           (sigma_to_med gt 1.) and (sigma_to_med lt 5) and (super_im_med_all lt 0.1) and (super_im_med_all gt 0.01): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (0.3*super_im_sig_all)), ocount_p)
           (sigma_to_med gt 1.) and (sigma_to_med lt 5) and (super_im_med_all lt 0.01): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (3*super_im_sig_all)), ocount_p)
           (sigma_to_med gt 5.): ptr_p = where(super_im_med_col_vec lt (super_im_med_all - (5*super_im_sig_all)), ocount_p)

        endcase 


               if (ocount_p ne 0) then begin 

                xp=ptr_p
                sptr_p = bsort(xp)
                xp = xp[sptr_p]
                sptr_p = uniq(xp)
                xp = xp[sptr_p]
                ;print,'xp =',xp
                print,''

                ymax_col = fltarr(ocount_p) 

                for z =0,ocount_p-1 do begin 

                 ymax_col[z] = max(image(xp[z],*),/NaN)

                endfor           

                ymax_col = ymax_col(where(finite(ymax_col)))
                pthres = min(ymax_col)/ft
                ;print,'pthres =',pthres

             endif else pthres = 0 

      endif

; Remove old flagging

        ptr = where(mask eq pulldownflag, ocount)
	if (ocount gt 0) then mask[ptr] = mask[ptr] - pulldownflag
	sxaddhist, name + ' Pulldown bit reset', hm


; If no pixels above pulldown threshold, inform user

       if (ocount_p eq 0) then begin
 
		mstr = 'No pulldown ' 
                mstr = mstr + ' for ' + strn(ifile) 
                print, mstr
                ;
                pthres = 0 & ft = 0 & backg = 0

        endif


  if (ocount_p ne 0) and (pthres gt 0) then begin

;   mask the columns affected by pulldown
	mask[xp, *] = mask[xp, *] or pulldownflag
	
; Next, we estimate what the image without artifacts should look like
; First generate boolean mask for pixels to interpolate
	bptr = where(mask and bad_pix_bits, bitcount)
	bmask = mask * 0B
	if (bitcount gt 0) then bmask[bptr] = 1B

; Weights for Gaussian smoothing kernel
        irac_p_make_gweights, 5, 5, sgauss[ch-1], w5
        irac_p_make_gweights, 11, 11, sgauss[ch-1], w11

; Once the flagging has been performed, fill in the masked pixels using
; some estimation scheme, usually a local interpolation but other schemes
; are possible too


; This is the default method using a bootstrapped Gaussian interpolate 
; and a last pass with a larger kernal
	irac_p_estimate_sky, image, bmask, w5, w11, estimated, 	/BOOTSTRAP, $
	                     FALLBACK_INTERP='NAN

; compute estimate of background in interpolated image 

        backg = median(estimated)

	
; Create difference image should just contain the artifacts, if the estimation 
; is good.
	diff = image - estimated

; Now fit the column pulldown, using the flagging from above, instead of
; re-flagging from the mask
		
; 
	for i = 0,ocount_p-1 do begin

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


; Find maximum pixel in column, this is the dominant pulldown trigger
		cptr = where(vec eq max(vec, /NAN))
; Locate row of maximum pixel in column		
		ymax = index[cptr[0]]
		
; Count top 5 pixels as a metric -- not used but could be if we try an apriori 
; solution, heritage code at this point
		sptr = bsort(vec[gptr], /REVERSE)
		top5 = total(vec[gptr[sptr[0:4]]])		
		
; Make separate corrections above and below trigger
; Top of column
		topptr = where(index gt ymax, tcount)
; Bottom of column
		botptr = where(index lt ymax, bcount)
		
; Initialize correction vector
		corr_vec = vec * 0.0	
; If enough pixel above source exist, then fit a 2nd order polynomial to them, 
; require polynomial order + 2 good data points

; Work on top first
               
               if (tcount ne 0) then begin

                tgptr = where(vec[topptr] eq vec[topptr], tgcount)
                if (tgcount gt npoints) then begin
; In addition should probably clip outliers before fitting: apply iterative sigma clipping 

                  ;  First, do basic statistics       
                  m0_top = median(vec[topptr[tgptr]]) 
                  mom_top = moment(vec[topptr[tgptr]], /NAN)
                  osig_top = sqrt(mom_top[1])                     
   
                  ; Then, do robust estimation of sigma       
                  sig_top = robust_sigma(vec[topptr[tgptr]])  
                  ; Uncertainty in median is estimated as uncertainty in mean  
                  sm0_top = sig_top / sqrt(tgcount)                            
                             
                  ; now identify and keep only values within a given sigma : iterative sigma-clipping

                   sig_thr_top = 30 
                   top_ct_pos = 0 & im_top_ct_neg = 0
                   tgptr_3s = where(vec[topptr[tgptr]] gt m0_top-sig_thr_top*sm0_top and vec[topptr[tgptr]] lt  m0_top+sig_thr_top*sm0_top, tp3scount)
                   mt = median(vec[topptr[tgptr[tgptr_3s]]])
                   pt0 = [mt, dblarr(npoints)]
                   ftt = {FLUX:vec[topptr[tgptr[tgptr_3s]]], X:index[topptr[tgptr[tgptr_3s]]]-ymax, $
                               ERR:svec[topptr[tgptr[tgptr_3s]]]}

                   pt = mpfit('exp_func_for_mpfit', pt0, FUNCTARGS=ftt, /quiet)
                   corr_vec[topptr[tgptr]] = pt(0) + pt(1) * exp(-(index[topptr[tgptr]]-ymax)^pt(2))
                   top_corr_pos =  where((corr_vec[topptr[tgptr]]) > 0, top_ct_pos)
                   im_top_corr = where(vec_im[topptr[tgptr]] - corr_vec[topptr[tgptr]] < 0, im_top_ct_neg)


                   cnt_top = 0
                   while (top_ct_pos eq 0) and (im_top_ct_neg eq 0)  and (sig_thr_top ge 3) do begin

                      sig_thr_top = sig_thr_top - cnt_top
                      tgptr_3s = where(vec[topptr[tgptr]] gt m0_top-sig_thr_top*sm0_top and vec[topptr[tgptr]] lt  m0_top+sig_thr_top*sm0_top, tp3scount)
                      mt = median(vec[topptr[tgptr[tgptr_3s]]])
                      pt0 = [mt, dblarr(npoints)]
                      ftt = {FLUX:vec[topptr[tgptr[tgptr_3s]]], X:index[topptr[tgptr[tgptr_3s]]]-ymax, $
                               ERR:svec[topptr[tgptr[tgptr_3s]]]}

                       pt = mpfit('exp_func_for_mpfit', pt0, FUNCTARGS=ftt, /quiet)
                       corr_vec[topptr[tgptr]] = pt(0) + pt(1) * exp(-(index[topptr[tgptr]]-ymax)^pt(2))
                       top_corr_pos =  where((corr_vec[topptr[tgptr]]) > 0, top_ct_pos)
                       im_top_corr = where(vec_im[topptr[tgptr]] - corr_vec[topptr[tgptr]] < 0, im_top_ct_neg)
                       cnt_top = cnt_top+1

                   endwhile

; If not enough good data then default to median of top portion for correction
                endif else begin   
                   mt = median(vec[topptr])   
                   if (tgcount gt 0) then begin 

                    corr_vec[tgptr] = mt
                    pt = fltarr(3)
                    pt(0) = mt
                    pt(1) = 0
                    pt(2) = 0

                   endif
 
               endelse
             
            endif

        
		
; Now same procedure on bottom of column	
 
                if (bcount ne 0) then begin 

                 bgptr = where(vec[botptr] eq vec[botptr], bgcount)
		 if (bgcount gt npoints) then begin
; In addition clip outliers before fitting: apply iterative sigma clipping 

                  ;  First, do basic statistics            
                  m0_bot = median(vec[botptr[bgptr]])      
                  mom_bot = moment(vec[botptr[bgptr]], /NAN)   
                  osig_bot = sqrt(mom_bot[1])                  
   
                  ; Then, do robust estimation of sigma        
                  sig_bot = robust_sigma(vec[botptr[bgptr]])    
                  ; Uncertainty in median is estimated as uncertainty in mean   
                  sm0_bot = sig_bot / sqrt(bgcount)                             

                  ; now keep only values > or < a given sigma : : iterative sigma-clipping

                   sig_thr_bot = 30
                   bot_ct_pos = 0 & im_bot_ct_neg = 0
                   bgptr_3s = where(vec[botptr[bgptr]] gt m0_bot-sig_thr_bot*sm0_bot and vec[botptr[bgptr]] lt  m0_bot + sig_thr_bot*sm0_bot, bt3scount)
                   mb = median(vec[botptr[bgptr[bgptr_3s]]])
                   pb0 = [mb, dblarr(npoints)] 
                   fb = {FLUX:vec[botptr[bgptr[bgptr_3s]]], X:ymax-index[botptr[bgptr[bgptr_3s]]], $
                               ERR:svec[botptr[bgptr[bgptr_3s]]]}

                   pb = mpfit('exp_func_for_mpfit', pb0, FUNCTARGS=fb,/quiet)
                   corr_vec[botptr[bgptr]] = pb(0) + pb(1) * exp(-(ymax-index[botptr[bgptr]])^pb(2))
                   bot_corr_pos =  where((corr_vec[botptr[bgptr]]) > 0, bot_ct_pos)
                   im_bot_corr = where(vec_im[botptr[bgptr]] - corr_vec[botptr[bgptr]] < 0, im_bot_ct_neg)




                   cnt_bot = 0
                   while (bot_ct_pos eq 0) and (im_bot_ct_neg eq 0)  and (sig_thr_bot ge 3) do begin

                      sig_thr_bot = sig_thr_bot - cnt_bot
                      bgptr_3s = where(vec[botptr[bgptr]] gt m0_bot-sig_thr_bot*sm0_bot and vec[botptr[bgptr]] lt  m0_bot + sig_thr_bot*sm0_bot, bt3scount)
                      mb = median(vec[botptr[bgptr[bgptr_3s]]])
                      pb0 = [mb, dblarr(npoints)]
                      fb = {FLUX:vec[botptr[bgptr[bgptr_3s]]], X:ymax-index[botptr[bgptr[bgptr_3s]]], $
                               ERR:svec[botptr[bgptr[bgptr_3s]]]}

                      pb = mpfit('exp_func_for_mpfit', pb0, FUNCTARGS=fb,/quiet)
                      corr_vec[botptr[bgptr]] = pb(0) + pb(1) * exp(-(ymax-index[botptr[bgptr]])^pb(2))
                      bot_corr_pos =  where((corr_vec[botptr[bgptr]]) > 0, bot_ct_pos)
                      im_bot_corr = where(vec_im[botptr[bgptr]] - corr_vec[botptr[bgptr]] < 0, im_bot_ct_neg)
                      cnt_bot = cnt_bot+1

                   endwhile


; If not enough good data then default to median of top portion for correction
		endif else begin 
            
                  mb = median(vec[botptr]) 
                  if (bgcount gt 0) then begin 

                   corr_vec[bgptr] = mb
                   pb = fltarr(3)
                   pb(0) = mb
                   pb(1) = 0
                   pb(2) = 0

                  endif

                endelse

              endif


		
; Now average adjacent rows' solutions to fit trigger pixel	
		if (ymax eq 0) then corr_vec[ymax] = corr_vec[ymax + 1] else $
		if (ymax eq ny-1) then corr_vec[ymax] = corr_vec[ymax - 1] else $ 	
		corr_vec[ymax] = (corr_vec[ymax + 1] + corr_vec[ymax - 1]) / 2.0
		
; Update image with correction - interpolate NaN values in fit-solution

                ptr_NaN = where(finite(corr_vec,/NaN), count_nan) 
                if (count_nan ne 0) then begin 
 
                 corr_vec_inter = interpol(corr_vec, findgen(n_elements(corr_vec)) , ptr_NaN)
                 corr_vec(ptr_NaN) = corr_vec_inter

                endif
		image[xp[i], *] = image[xp[i], *] - corr_vec
; and update sigma
		sigma[xp[i], *] = sigma[xp[i], *] * pullnoise_factor
		


;          endif

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
	
   

 endif else begin 

        writefits, ofile, image, h
        writefits, osfile, sigma, hs
        writefits, omfile, mask, hm

 endelse

; endif


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
                        forward_function poly2d_func 
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


PRO MEDCLIP, Image, Mean, Sigma, CLIPSIG=clipsig, MAXITER=maxiter, $
    CONVERGE_NUM=converge_num, VERBOSE=verbose, SUBS=subs
;+
; NAME:
;       MEDCLIP
;
; PURPOSE:
;       Computes an iteratively sigma-clipped mean on a data set
; EXPLANATION:
;       Clipping is done about median, but mean is returned.
;       Called by SKYADJ_CUBE
;
; CATEGORY:
;       Statistics
;
; CALLING SEQUENCE:
;       MEANCLIP, Data, Mean, Sigma
;
; INPUT POSITIONAL PARAMETERS:
;       Data:     Input data, any numeric array
;
; OUTPUT POSITIONAL PARAMETERS:
;       Mean:     N-sigma clipped mean.
;       Sigma:    Standard deviation of remaining pixels.
;
; INPUT KEYWORD PARAMETERS:
;       CLIPSIG:  Number of sigma at which to clip.  Default=3
;       MAXITER:  Ceiling on number of clipping iterations.  Default=5
;       CONVERGE_NUM:  If the proportion of rejected pixels is less
;           than this fraction, the iterations stop.  Default=0.02, i.e.,
;           iteration stops if fewer than 2% of pixels excluded.
;       /VERBOSE:  Set this flag to get messages.
;
; OUTPUT KEYWORD PARAMETER:
;       SUBS:     Subscript array for pixels finally used.
;
;
; MODIFICATION HISTORY:
;       Based on meanclip written by:   RSH, RITSS, 21 Oct 98
;       20 Jan 99 - Added SUBS, fixed misplaced paren on float call,
;                   improved doc.  RSH
;       Edited by R. Paladini, 14 Aug 2010 to output median instead of mean
;-

IF n_params(0) LT 1 THEN BEGIN
    print, 'CALLING SEQUENCE:  MEANCLIP, Image, Mean, Sigma'
    print, 'KEYWORD PARAMETERS:  CLIPSIG, MAXITER, CONVERGE_NUM, ' $
        + 'VERBOSE, SUBS'
    RETALL
ENDIF

prf = 'MEANCLIP:  '

verbose = keyword_set(verbose)
IF n_elements(maxiter) LT 1 THEN maxiter = 5
IF n_elements(clipsig) LT 1 THEN clipsig = 3
IF n_elements(converge_num) LT 1 THEN converge_num = 0.02

subs = where(finite(image),ct)
iter=0
REPEAT BEGIN
    skpix = image[subs]
    iter = iter + 1
    lastct = ct
    medval = median(skpix)
    sig = stdev(skpix)
    wsm = where(abs(skpix-medval) LT clipsig*sig,ct)
    IF ct GT 0 THEN subs = subs[wsm]
ENDREP UNTIL (float(abs(ct-lastct))/lastct LT converge_num) $
          OR (iter GT maxiter)


med = median(image[subs])
mean = (moment(med))[0]
sigma = stddev(image[subs],/NaN)


IF verbose THEN BEGIN
    print, prf+strn(clipsig)+'-sigma clipped mean'
    print, prf+'Mean computed in ',iter,' iterations'
    print, prf+'Mean = ',mean,',  sigma = ',sigma
ENDIF

RETURN
END

