pro mips24_lir

readcol, '/Users/jkrick/nep/clusters/mips24/bgs.txt', number,  f12,  f25,  f60,  f100,  cz, dist, loglir1, loglir2, format='F,F,F,F,F,F,F,F,F',/debug,skipline=2

;lir1 = total far infrared luminisuty 40-400micron
;lir2 = total infrared luminosity 8-1000micron

lir1 = 10^(loglir1)
lir2 = 10^(loglir2)

!p.multi=[0,1,2]
plot, f12, loglir1, psym = 2,xrange=[1E-2,200],/xlog,/ylog
plot, f12, loglir2, psym = 3, xrange=[1E-2,100],/xlog,yrange=[8,13]

;fit  the relation for lir2
start = [-2.5, 12] 
result = MPFITFUN('linear', f12, loglir2, err, start)

oplot, findgen(100)/10, result(0)*findgen(100)/10 + result(1)
oplot, findgen(100)/10, -2.5*findgen(100)/10 + 12
;oplot, f12, -2.5*f12 + 15
;oplot, [.2, 1.0], [10^(11.5), 1E10], thick = 3


oplot, findgen(100)/10, 0.89*(findgen(100)/10)^(1.094)

end
