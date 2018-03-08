function obsv_filter,inlam,innulnu,$
mips24=mips24

; function that will calculate the observed lnu in a given filter
; window. Input is lambda (in microns), nulnu in say Lsun; output is lnu in Lsun/Hz
; Assumes filter curves are in microns and r_lambda

COMMON FILTERCURVE

if keyword_set(mips24) then  begin
  flb = filters[27].lamb 
  transp = filters[27].transp
  goto,jumpt
endif

jumpt: if keyword_set(boxcar) then begin
   flb = lrange[0] + findgen(101)*(lrange[1]-lrange[0])/100.
   transp = fltarr(101)+1.
endif else begin
   lrange = minmax(flb)
endelse

delta_lam = (lrange[1]-lrange[0])/100.
goodrng = where(inlam ge lrange[0] and inlam le lrange[1])

if (goodrng[0] ne -1) then begin
  ltz = where(innulnu lt 0.)
  if (ltz[0] ne -1) then begin
    ;print,'WARNING: SED has negative values.'
    innulnu[ltz] = 0.
  endif
  inlnu = innulnu/(3d14/inlam)

  newinlnu = interpol(inlnu,inlam,flb)
  notval = where((flb lt min(inlam)) or (flb gt max(inlam)))
  if (notval[0] ne -1) then newinlnu[notval] = 0.

  dfreq = 3d14/(flb-0.5*delta_lam)-3d14/(flb+0.5*delta_lam)
 
  totval = total(newinlnu*transp*dfreq)/total(transp*dfreq)

  return,totval[0]
endif else return,0.     ; if filter doesnt see the flux, it is 0

end
