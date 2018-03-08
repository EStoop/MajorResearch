pro avg_stack
  restore, '/Users/jkrick/irac_warm/snapshots/hd7924/hd7924_old_corr.sav'
  vec = snapshots[0].flux
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

;why are there negative values?
  good = where(db gt 0)
  db = db[good]
  bs = bs[good]
  nwins = n_elements(db)
;change variable names to what makes sense to me:
  print, 'nwins', nwins
  print, 'db', db
  print, 'bs', bs
  nfits = db 
  nframes = bs

  binnedflux = fltarr(nwins)
  binnedfluxerr = fltarr(nwins)
  binnedcorrflux = fltarr(nwins)

  for nw = 0, nwins - 1 do begin  ; for each binning size
     ;print, 'nw', nw
     for si = 0., nfits[nw] - 1 do begin  ; for the number of elements of that bin
        ;uncorrected flux
        idata = snapshots[0].flux[si*nframes[nw]:si*nframes[nw] + (nframes[nw] - 1)]
        idataerr = snapshots[0].fluxerr[si*nframes[nw]:si*nframes[nw] + (nframes[nw] - 1)]
        binnedflux[nw] = mean(idata,/nan)
        binnedfluxerr[nw] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        
        ;corrected flux
        cdata = snapshots[0].corrflux[si*nframes[nw]:si*nframes[nw] + (nframes[nw] - 1)]
        binnedcorrflux[nw] = mean(cdata,/nan)
 
     endfor
  endfor
  
  binnedflux = binnedflux ;/ binnedflux(0)
  binnedcorrflux = binnedcorrflux ;/ binnedcorrflux(0)
  
  a = plot( nfits, binnedflux,/xlog,/ylog, '1o', sym_size = 0.6, sym_filled = 1, color='black', name = 'original flux');, yrange = [0.1, 1.0])
  b = plot(nfits, binnedcorrflux, '1s', color = 'red',/overplot, name = 'pmap corrected')
  l =   legend(Target = [a,b], position = [1E3,0.8], /data)
  sroot = sqrt(nfits)
  c = plot(nfits, 1/(sroot),/overplot)
  
end
