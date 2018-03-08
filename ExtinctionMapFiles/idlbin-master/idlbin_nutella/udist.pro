pro udist

restore, '/Users/jkrick/idlbin/object.sav'

ps_open, filename='/Users/jkrick/palomar/lfc/catalog/udist.ps',/portrait,/square,/color

plothist, object.umaga, xhist, yhist, bin=0.1, /noprint, xrange=[15,30], thick = 3, xthick=3, ythick = 3, charthick = 3, xtitle="u' AB magnitude", ytitle = 'N'

ps_close, /noprint, /noid
end
