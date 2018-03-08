pro avg_stack, vec, db, rs, ds, mw, sw
; GOO need to revise for current purposes!!!
	nvec = n_elements(vec)

; Average progressively larger windows, calculate mean and dispersion for each one,
; find all windows such that n is modulo window
	n4 = nvec / 4
; Start with all possible bins
	db = indgen(n4) + 2
	bs = nvec / db
; Find bin sizes that are evenly divisible	
	nrem = nvec mod db
	ptr = where(nrem eq 0, count)
	if (count eq 0) then message, 'No windows are even divisors of sample'
; Keep the evenly divisible bins
	bs = bs[ptr]
	db = db[ptr]
	nwins = n_elements(bs)
	mw = dblarr(nwins)
	sw = dblarr(nwins)
	rs = sw
	ds = sw
	
;	bvec = dblarr(bs[0], nwins) + !VALUES.D_NAN

; Loop over all window sizes, smoothing data and calculating the mean and dispersion.

	for k = 0, nwins-1 do begin
		tvec = dblarr(bs[k])
		for l = 0, bs[k]-1 do tvec[l] = total(vec[l*db[k]:(l+1)*db[k]-1]) / db[k]
;		bvec[0:bs[k]-1, k] = tvec
		mom = moment(tvec)
		mw[k] = mom[0]
		sw[k] = sqrt(mom[1])
	endfor
	rs = sw / sw[0]
	ds = rs - 1. / sqrt(db / db[0])

return
end