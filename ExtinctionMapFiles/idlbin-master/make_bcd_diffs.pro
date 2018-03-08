pro make_bcd_diffs, aor, channel, ADJUST=adjust, WRITE_RATIO=write_ratio
	if (N_params() ne 2) then begin
		print, 'Syntax == make_bcd_diffs, reqkey, channel, /ADJUST, /WRITE_RATIO'
		return
	endif
	

	if (not keyword_set(ADJUST)) then adjust=0
	thresh = 3.0D-05
; Code if you need or want to remove flatfields
	get_flats, flata, flatb, sflata, sflatb
	mflatrat = fltarr(4)
	msflatrat = fltarr(4)
	for i = 0, 3 do begin
		mflatrat[i] = mean(flata[*, *, i] / flatb[*, *, i])
		msflatrat[i] = mean(sflata[*, *, i] / sflatb[*, *, i])
	endfor
;	print, mflatrat
;	print, msflatrat
	if (n_elements(flata) eq 1) then adjust=0
	
	dir = '~/iracdata/FINAL_CHECK/IRAC012200/'
	ndir = dir + 'S18.18v2/'
	odir = dir + 'S18.7/'
	
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
	xstr = odir + suffix
	print, 'Processing request key ' + saor

        ps_start,filename = sstr +  '/diff_histograms.ps'


	files = file_search(sstr + '/SPITZER_I' + sch + '_' + saor + '_*_bcd.fits', COUNT=nfiles)
	xfiles = file_search(xstr + '/SPITZER_I' + sch + '_' + saor + '_*_bcd.fits', COUNT=xcount)
	if (xcount eq nfiles) then begin
		for j = 0, nfiles-1 do begin
			print_pct, j, nfiles

			nim = readfits(files[j], nh, /SILENT)
			flist[k] = files[j]
			remove_leading_path, files[j], file
			fsrch = xstr + '/' + strmid(file, 0, 30) + '*_bcd.fits'
			ofiles = file_search(fsrch, COUNT=ocount)
			if (ocount eq 1) then ok++ else begin
				if (ocount eq 0) then print, 'No corresponding file for ' + file
				if (ocount gt 0) then print, 'Multiple files for ' + file
				stop
			endelse
			oim = readfits(ofiles[0], oh, /SILENT)
			
			nch = sxpar(nh, 'CHNLNUM')
			ch[k] = nch
			nfc = sxpar(nh, 'FLUXCONV')
			ofc = sxpar(oh, 'FLUXCONV')
			is_sub = sxpar(nh, 'AREADMOD')		
			
			if (adjust) then begin
				xim = oim * nfc / ofc 
				if (is_sub) then for sf = 0, 63 do $
						xim[*, *, sf] = xim[*, *, sf] * sflatb[*, *,nch-1] / sflata[*,*,nch-1] $ 
				else xim = xim * flatb[*, *, nch-1] / flata[*, *, nch-1]
			endif else xim = oim
			diff = xim - nim

			dh = nh
			sxaddpar, dh, 'BITPIX', -64, ' Double precision'
			dmin = min(diff,/NAN)
			dmax = max(diff,/NAN)
			dmean = mean(diff, /NAN)
                        
                        plothist, diff, xhist, yhist, /nan, bin = 0.00005, xrange = [-.005,0.005],/ylog, xtitle = 'diff value', ytitle = 'Number', title = 'frame' + string(j)

			if (adjust) then begin
				sxaddpar, dh, 'OFLUXCNV', ofc, ' Old flux conversion'
				if (is_sub) then sxaddpar, dh, 'MFLATRAT', msflatrat[nch-1], ' Flat avg. diff' $
							else sxaddpar, dh, 'MFLATRAT', mflatrat[nch-1], ' Flat avg. diff'
				sxaddpar, dh, 'DIFFTYPE', 'ADJUST', ' flat-field, fluxconv backed out'		
			endif else begin
				sxaddpar, dh, 'DIFFTYPE', 'RAW', ' raw difference, does not backout flat'
			endelse
			sxaddpar, dh, 'FTYPE', 'diff', ' Difference between old and new processing'
			sxaddpar, dh, 'DMIN', dmin, ' Minimum difference between old and new'
			sxaddpar, dh, 'DMAX', dmax, ' Maximum difference between old and new'
			sxaddpar, dh, 'DMEAN', dmean, ' Mean difference between old and new'

			rh = dh
			rat = diff / nim
			ptr = where(abs(rat) gt thresh, tcount)
			sxaddpar, rh, 'THRES', thresh, ' Threshold for outliers in ratio'
			sxaddpar, rh, 'OCOUNT', tcount, ' Outliers in ratio'
			ptr = where(abs(rat) gt 0.01, tcount1)
			sxaddpar, rh, 'OCOUNT', tcount1, ' Outliers greater than 0.1%'
			
			rename, files[j], '_bcd', '_diff', dfile
			writefits, dfile, diff, dh
			rename, files[j], '_bcd', '_rat', rfile
			if (keyword_set(WRITE_RATIO)) then writefits, rfile, rat, rh	
			
			nthresh[k] = tcount
			nthresh1[k] = tcount1
			mdiff[k] = dmean
			k = k + 1L
		endfor
	endif
	
	nthresh = nthresh[0L:k-1L]
	nthresh1 = nthresh1[0L:k-1L]
	flist = flist[0L:k-1L]
	mdiff = mdiff[0L:k-1L]	
	ch = ch[0L:k-1L]
	print, 'Number of new files = ' + strn(k)
	print, 'Number of old files = ' + strn(ok)
	save, FILENAME=sstr+'bcd_diff.sav', nthresh, nthresh1, mdiff, flist, ch

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
