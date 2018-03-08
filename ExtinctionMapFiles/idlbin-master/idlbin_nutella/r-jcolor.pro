pro rjcolor

restore, '/Users/jkrick/idlbin/object.sav'
good = where(object.rmaga gt 0 and object.rmaga lt 90 and object.wircjmag gt 0 and object.wircjmag lt 90 and object.zphot lt 1.2 and object.zphot gt 0.8 and object.spt lt 3)


wircjab = 1594.*10^(object[good].wircjmag/(-2.5))
wircjmagab = -2.5*alog10(wircjab) +8.926

plothist, object[good].rmaga - wircjmagab, thick=1, bin=0.1, charthick=1


end
