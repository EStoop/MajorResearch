pro sed

close, /all
device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

readcol,'/Users/jkrick/ZPHOT/templates/e0.1.txt',lambda1, j,j,j,j,flux1,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e0.2.txt',lambda2, j,j,j,j,flux2,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e0.4.txt',lambda4, j,j,j,j,flux4,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e0.8.txt',lambda8, j,j,j,j,flux8,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e1.5.txt',lambda15, j,j,j,j,flux15,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e2.0.txt',lambda20, j,j,j,j,flux20,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e3.0.txt',lambda30, j,j,j,j,flux30,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e4.0.txt',lambda40, j,j,j,j,flux40,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e5.0.txt',lambda50, j,j,j,j,flux50,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e8.0.txt',lambda80, j,j,j,j,flux80,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e11.0.txt',lambda110, j,j,j,j,flux110,j,j,j,format="A"
readcol,'/Users/jkrick/ZPHOT/templates/e13.0.txt',lambda130, j,j,j,j,flux130,j,j,j,format="A"


ps_open, filename='/Users/jkrick/palomar/lfc/catalog/ellip.sed.ps',/portrait,/square,/color

plot, lambda1, flux1, thick = 3, xthick=3, ythick=3,/xlog,/ylog, color = colors.black
oplot, lambda2, flux2, thick = 3, color = colors.black
oplot, lambda4, flux4, thick = 3, color = colors.gray
oplot, lambda8, flux8, thick = 3, color = colors.brown
oplot, lambda15, flux15, thick = 3, color = colors.red
oplot, lambda20, flux20, thick = 3, color = colors.orange
oplot, lambda30, flux30, thick = 3, color = colors.green
oplot, lambda40, flux40, thick = 3, color = colors.cyan
oplot, lambda50, flux50, thick = 3, color = colors.blue
oplot, lambda80, flux80, thick = 3, color = colors.purple
oplot, lambda110, flux110, thick = 3, color = colors.violet
oplot, lambda130, flux130, thick = 3, color = colors.pink


ps_close, /noprint,/noid

end
