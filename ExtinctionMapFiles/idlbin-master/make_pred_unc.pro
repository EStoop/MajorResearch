pro make_pred_unc, aor, channel, WRITE_RATIO=write_ratio
;+
; NAME: 
;      MAKE_PRED_UNC
;
; PURPOSE:
;      To make an estimate of the BCD uncertainty from the input BCD.  Uncertainty is
;      assumed to be (Poisson^2 + readnoise^2)^0.5 where the Poisson noise is estimated from
;      the BCD directly.
;
; HISTORY:
;      Fixed conversion back to MJy/sr 14 Oct 2010 SJC
;-

; Modify path to iracdata if needed
	maindir = '~/iracdata/'
	
	if (N_params() ne 2) then begin
		print, 'Syntax == make_pred_unc, reqkey, channel, /WRITE_RATIO'
		return
	endif
	
	thresh = 3.0D-05	
; Code if you need or want to remove flatfields
	get_flats, flata, flatb, sflata, sflatb	
		
	dir = maindir + 'FINAL_CHECK/IRAC012200/'
	ndir = dir + 'S18.18v2/'

; Readnoise for each frametime, each channel some frametimes are not used in cryogenic mission
	readnoise = [[23.9, 23.8, 12.7, 11.9], [16.9, 16.8, 9.0, 8.4], $
			 [11.8, 12.1, 9.1, 7.1], [9.4, 9.4, 8.8, 6.7], $
			 [23.9, 23.8, 12.7, 11.9], [23.9, 23.8, 12.7, 11.9], $
			 [23.9, 23.8, 12.7, 11.9], [11.8, 12.1, 9.1, 7.1], $
			 [9.4, 9.4, 8.8, 6.7], [9.4, 9.4, 8.8, 6.7], $
			 [7.8, 7.5, 10.7, 6.9], [8.4, 7.9, 13.1, 6.8]] * 1.D	
	
	nthresh = fltarr(200000L)
	nthresh1 = fltarr(200000L)
	mdiff = fltarr(200000L)
	flist = strarr(200000L)
	ch = fltarr(200000L)
	
	k = 0L
	ok = 0L
	
	saor = strn(aor)
	sch = strn(channel)
	suffix = 'r' + saor + '/ch' + sch + '/bcd/'
	sstr = ndir + suffix
	print, 'Processing request key ' + saor

	files = file_search(sstr + '/SPITZER_I' + sch + '_' + saor + '_*_bcd.fits', COUNT=nfiles)
        ps_start,filename = sstr +  '/unc_histograms.ps'

	for j = 0, nfiles-1 do begin
		print_pct, j, nfiles

		nim = readfits(files[j], nh, /SILENT)
		rename, files[j], '_bcd', '_bunc', ufile
		unc = readfits(ufile, uh, /SILENT)
		flist[k] = files[j]
		remove_leading_path, files[j], file
		
;		rename, istr, oldsub, newsub, nstr

		nch = sxpar(nh, 'CHNLNUM')
		ch[k] = nch
		fluxconv = sxpar(nh, 'FLUXCONV')
		gain = sxpar(nh, 'GAIN')
		areadmod = sxpar(nh, 'AREADMOD')	
		ft = sxpar(nh, 'FRAMTIME')
		tex = sxpar(nh, 'EXPTIME')
		skydk = sxpar(nh, 'SKYDKMED')

; Get readnoise for observation
		if (abs(ft-0.02) lt 0.01) then findex = 0 $
		else if (abs(ft-0.1) lt 0.01) then findex = 1 $
		else if (abs(ft-0.4) lt 0.01 and areadmod eq 1) then findex = 2 $
		else if (abs(ft-0.4) lt 0.01 and areadmod eq 0) then findex = 4 $
		else if (abs(ft-0.6) lt 0.01) then findex = 5 $
		else if (abs(ft-1.2) lt 0.01) then findex = 6 $
		else if (abs(ft-2.) lt 0.01 and areadmod eq 1) then findex = 3 $
		else if (abs(ft-2.) lt 0.01 and areadmod eq 0) then findex = 7 $
		else if (abs(ft-6.) lt 0.01) then findex = 8 $
		else if (abs(ft-12.) lt 0.01) then findex = 9 $
		else if (abs(ft-30.) lt 0.01) then findex = 10 $
		else if (abs(ft-100.) lt 0.01) then findex = 11 else findex = 0		
; Readnoise for current BCD
		rn = readnoise[nch-1, findex]

; Convert image to equivalent electrons
; First remove flat
		nunc = nim * flata[*, *, nch-1]
; and add skydark
		nunc = nunc + skydk
; Nelectrons = MJy/sr * GAIN * EXPTIME / FLUXCONV 
		nunc = nunc * gain * tex / fluxconv
; add readnoise in quadrature 
		nunc = nunc + rn * rn
; return sqrt(variance)
		nunc = sqrt(nunc)
; Convert back to MJy/sr
		nunc = nunc * fluxconv / (gain * tex)
		
; Make new header		
		nuh = nh
		sxaddpar, nuh, 'BITPIX', -64, ' Double precision'
		sxaddpar, nuh, 'FTYPE', 'PREDUNC', ' Predicted uncertainty'		
		rename, files[j], '_bcd', '_punc', nufile
; Write predicted uncertainty
		writefits, nufile, nunc, nuh
		
; Make differnce between rough uncertainty and pipeline uncertainty	
		diff = unc - nunc

		dh = uh
		dmin = min(diff)
		dmax = max(diff)
		dmean = mean(diff, /NAN)
                plothist, diff, xhist, yhist, /nan, bin = 0.05, xrange = [-2,2],/ylog, xtitle = 'diff value', ytitle = 'Number', title = 'frame' + string(j)

		sxaddpar, dh, 'BITPIX', -64, ' Double precision'
		sxaddpar, dh, 'DIFFTYPE', 'RAW', ' raw difference, does not backout flat'
		sxaddpar, dh, 'FTYPE', 'uncdiff', ' Difference between pipeline unc and quick unc'
		sxaddpar, dh, 'DMIN', dmin, ' Minimum difference between uncertainties'
		sxaddpar, dh, 'DMAX', dmax, ' Maximum difference between uncertainties'
		sxaddpar, dh, 'DMEAN', dmean, ' Mean difference between uncertainties'
		rename, files[j], '_bcd', '_uncdiff', dfile		
		writefits, dfile, diff, dh

		rh = dh
		rat = diff / unc
		ptr = where(abs(rat) gt thresh, tcount)
		
		sxaddpar, rh, 'THRES', thresh, ' Threshold for outliers in ratio'
		sxaddpar, rh, 'OCOUNT', tcount, ' Outliers in ratio'
		ptr = where(abs(rat) gt 0.01, tcount1)
		sxaddpar, rh, 'OCOUNT', tcount1, ' Outliers greater than 0.1%'
		rename, files[j], '_bcd', '_uncrat', rfile
		if (keyword_set(WRITE_RATIO)) then writefits, rfile, rat, rh
		
		nthresh[k] = tcount
		nthresh1[k] = tcount1
		mdiff[k] = dmean
		k = k + 1L

	endfor

; statistics for save file	
	nthresh = nthresh[0L:k-1L]
	nthresh1 = nthresh1[0L:k-1L]
	flist = flist[0L:k-1L]
	mdiff = mdiff[0L:k-1L]	
	ch = ch[0L:k-1L]
	save, FILENAME=sstr+'pred_unc.sav', nthresh, nthresh1, mdiff, flist, ch
        ps_end
	
return
end

;;pro rename, istr, oldsub, newsub, nstr
;;pro remove_leading_path, istr, ostr, LPATH=lpath
pro get_flats, flata, flatb, sflata, sflatb
	aflatdir = '~/iracdata/carey/cryo_abs_cal/new_cals/'
	
	flat1a = readfits(aflatdir + 'irac_b1_fa_superskyflat_finalcryo_021110.fits', hf1a, /SILENT)
	flat1a = flat1a[*, *, 0]
	flat1a = rotate(flat1a, 7)
	flat2a = readfits(aflatdir + 'irac_b2_fa_superskyflat_finalcryo_021110.fits', hf2a, /SILENT)
	flat2a = flat2a[*, *, 0]
	flat2a = rotate(flat2a, 7)
	flat3a = readfits(aflatdir + 'irac_b3_fa_superskyflat_finalcryo_021110.fits', hf3a, /SILENT)
	flat3a = flat3a[*, *, 0]
	flat4a = readfits(aflatdir + 'irac_b4_fa_superskyflat_finalcryo_021110.fits', hf4a, /SILENT)
	flat4a = flat4a[*, *, 0]
	
	sflat1a = readfits(aflatdir + 'irac_b1_sa_superskyflat_finalcryo_021110.fits', hsf1a, /SILENT)
	sflat1a = sflat1a[*, *, 0]
	sflat1a = rotate(sflat1a, 7)
	sflat2a = readfits(aflatdir + 'irac_b2_sa_superskyflat_finalcryo_021110.fits', hsf2a, /SILENT)
	sflat2a = sflat2a[*, *, 0]
	sflat2a = rotate(sflat2a, 7)
	sflat3a = readfits(aflatdir + 'irac_b3_sa_superskyflat_finalcryo_021110.fits', hsf3a, /SILENT)
	sflat3a = sflat3a[*, *, 0]
	sflat4a = readfits(aflatdir + 'irac_b4_sa_superskyflat_finalcryo_021110.fits', hsf4a, /SILENT)
	sflat4a = sflat4a[*, *, 0]
	
	bflatdir = '~/iracdata/carey/cryo_abs_cal/old_cals/'
	
	flat1b = readfits(bflatdir + 'irac_b1_fa_superskyflat_121407.fits', hfb1, /SILENT)
	flat1b = flat1b[*, *, 0]
	flat1b = rotate(flat1b, 7)
	flat2b = readfits(bflatdir + 'irac_b2_fa_superskyflat_121407.fits', hfb2, /SILENT)
	flat2b = flat2b[*, *, 0]
	flat2b = rotate(flat2b, 7)
	flat3b = readfits(bflatdir + 'irac_b3_fa_superskyflat_121407.fits', hfb3, /SILENT)
	flat3b = flat3b[*, *, 0]
	flat4b = readfits(bflatdir + 'irac_b4_fa_superskyflat_121407.fits', hfb4, /SILENT)
	flat4b = flat4b[*, *, 0]
	
	sflat1b = readfits(bflatdir + 'irac_b1_sa_superskyflat_121407.fits', hsfb1, /SILENT)
	sflat1b = sflat1b[*, *, 0]
	sflat1b = rotate(sflat1b, 7)
	sflat2b = readfits(bflatdir + 'irac_b2_sa_superskyflat_121407.fits', hsfb2, /SILENT)
	sflat2b = sflat2b[*, *, 0]
	sflat2b = rotate(sflat2b, 7)
	sflat3b = readfits(bflatdir + 'irac_b3_sa_superskyflat_121407.fits', hsfb3, /SILENT)
	sflat3b = sflat3b[*, *, 0]
	sflat4b = readfits(bflatdir + 'irac_b4_sa_superskyflat_121407.fits', hsfb4, /SILENT)
	sflat4b = sflat4b[*, *, 0]	
	
	flata = dblarr(256, 256, 4)
	flatb = flata
	sflata = dblarr(32, 32, 4)
	sflatb = dblarr(32, 32, 4)

	flata[*, *, 0] = flat1a
	flata[*, *, 1] = flat2a
	flata[*, *, 2] = flat3a
	flata[*, *, 3] = flat4a	
	
	flatb[*, *, 0] = flat1b
	flatb[*, *, 1] = flat2b
	flatb[*, *, 2] = flat3b
	flatb[*, *, 3] = flat4b	
	
	sflata[*, *, 0] = sflat1a
	sflata[*, *, 1] = sflat2a
	sflata[*, *, 2] = sflat3a
	sflata[*, *, 3] = sflat4a	
	
	sflatb[*, *, 0] = sflat1b
	sflatb[*, *, 1] = sflat2b
	sflatb[*, *, 2] = sflat3b
	sflatb[*, *, 3] = sflat4b	

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

pro rename, istr, oldsub, newsub, nstr
;+
; NAME:
;    RENAME
; 
; PURPOSE:
;    To change a filename by replacing a substring with a new substring
;
; INPUT:
;    istr: string containing input filename
;    oldsub: string containing substring to change
;    newsub: string containing substring to change to
;
; OUTPUT: 
;    nstr: string containing output filename
;    
; HISTORY
;    Initial coding SJC 15 Nov 2006
;    Module checked SJC 27 Dec 2006
;-
    nosub = strlen(oldsub)
    nstr = istr

    pos = strpos(nstr, oldsub)
    while (pos gt -1) do begin
        nnstr = strlen(nstr)
        nstr = strmid(nstr, 0, pos) + newsub + $
                        strmid(nstr, pos+nosub, nnstr-pos+nosub)
        pos = strpos(nstr, oldsub)
    endwhile

return
end

pro print_pct, count, max, pcts=pcts
	if (not keyword_set(PCTS)) then pcts = 10
	indices = (indgen(fix(max/pcts)>1)+1) * fix(pcts)
	vals = long(float(indices) / 100.0 * max)
	ptr = where(count eq vals, n)
	if (n eq 1) then print,FORMAT='($,"...",I3,"%")', indices[ptr[0]]
	if (count eq max-1) then begin 
		print, FORMAT='($,"...100%")'
		print, ''
	endif

return
end
