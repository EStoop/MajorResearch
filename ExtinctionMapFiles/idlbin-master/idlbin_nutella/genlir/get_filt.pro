pro get_filt

; function that will load all filter curves into memory
; Assumes filter curves are in microns and r_lambda

rootdir = './'
COMMON FILTERCURVE
filters = replicate({fc, lamb:fltarr(101), transp:fltarr(101)}, 32)

for i = 27,27 do begin

  case i of 
  27: readcol,'/Users/jkrick/nutella/idlbin/genlir/mips24lg.dat',flb,transp,format='D,D',/silent
  endcase

  maxtransp = max(transp)
  goodvals = where(transp ge 0.1*maxtransp)
  tlam = flb[goodvals]
  ttransp = transp[goodvals]
  lrange = minmax(tlam)

  flb = lrange[0] + findgen(101)*(lrange[1]-lrange[0])/100.
  transp = interpol(ttransp,tlam,flb)

  filters[i].lamb = flb
  filters[i].transp = transp

endfor

end
