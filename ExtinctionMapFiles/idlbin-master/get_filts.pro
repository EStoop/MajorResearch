pro get_filts

; function that will load all filter curves into memory
; Assumes filter curves are in microns and r_lambda

;rootdir = './'
rootdir = '~/Virgo/irac/'
COMMON FILTERCURVE, filters
filters = replicate({fc, lamb:fltarr(101), transp:fltarr(101)}, 2)

for i = 0,1 do begin

   case i of 
      0: readcol,rootdir+'ch1trans.txt',flb,transp,format='D,D',/silent
      1: readcol,rootdir+'ch2trans.txt',flb,transp,format='D,D',/silent
    endcase

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
