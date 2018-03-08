pro mode_calc, x, mode, VERBOSE=verbose
;+
; NAME:
;      MODE_CALC
;
; PURPOSE:
;      To calculate the mode of a distribution using the method of Bickel (2002)
;
; INPUT:
;      X: Array containing values to calculate mode of
;
; OUTPUT:
;      MODE: Mode of X
;
; FUNCTIONS USED:
;      MIN(), MAX(), FINITE(), BSORT()
;
; PROCEDURES USED:
;      PRINT_PCT
;
; HISTORY:
;      Header added, VERBOSE keyword added 25 Oct 2011 SJC
;
;-
; ALgorithm from Bickel (2002)

; Initialize data vector, remove NaNs, and sort
  v = x
  ptr = where(finite(v))
  v = v[ptr]
  n = n_elements(v)
  ptr = bsort(v)
  v = v[ptr]
  nmax = n
  nprime = 0L
  done = 0	
  while (not done) do begin
     if (keyword_set(VERBOSE)) then print, n, nmax
; If two or fewer samples left in vector, then average to get mode
     if (n lt 16) then begin
        mode = mean(v)
        done = 1
     endif
     
; Otherwise find range and...
     vmax = max(v)
     vmin = min(v)
     hrange = (vmax - vmin) / 2.
; Find number of samples in intervals v(i),v(i)+ hrange		
     ninterval = lonarr(n)
     interval_range = lonarr(n)
     for i = 0, n-1 do begin
        ptr = where(v ge v[i] and v le v[i] + hrange, count)
        ninterval[i] = count
        interval_range[i] = v[ptr[count-1]] - v[ptr[0]]
     endfor
     ptr = where(ninterval eq max(ninterval), count)
     
     if (count gt 1) then begin
        mrange = min(interval_range[ptr])
        ptr = where(ninterval eq max(ninterval) and $
                    interval_range eq mrange, count)
        vmax = max(v[ptr]) + hrange
        vmin = min(v[ptr])
        ptr = where(v ge vmin and v le vmax, nprime)
        if (keyword_set(VERBOSE)) then print, 'nprime', nprime
        
        v = v[ptr]
        if (nprime eq n and nprime gt 2) then begin
           vbot_diff = v[1] - v[0]
           vtop_diff = v[nprime-1L] - v[nprime-2L]
           if (vbot_diff lt vtop_diff) then v = v[0:nprime-2L] $
           else if (vbot_diff gt vtop_diff) then v = v[1:nprime-1L] $
           else v = v[1:nprime-2L]
        endif
        if (keyword_set(VERBOSE)) then print, 'nprime after', nprime

        n = n_elements(v)
     endif else begin
; If only one interval with max number of points, then restrict vector to that interval		
        vptr = where(v ge v[ptr[0]] and v le v[ptr[0]] + hrange, n)
        v = v[vptr]
     endelse		
  endwhile
  return
end
