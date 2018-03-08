pro phottest
colors = GetColor(/load, Start=1)

!p.multi = [0, 2, 1]
ps_open, file = "/Users/jkrick/palomar/lfc/sdss/phottest.ps", /portrait, xsize = 8, ysize = 4,/color

readcol, "/Users/jkrick/palomar/LFC/catalog/SExtractor.r.cat", o, xwcsr, ywcsr, xcenterr, ycenterr, fluxautor, magautor, magerrr, fwhmr, isoarear, fluxmaxr, ellipr
readcol, "/Users/jkrick/palomar/LFC/catalog/SExtractor.g.cat",o, xwcsg, ywcsg, xcenterg, ycenterg, fluxautog, magautog, magerrg, fwhmg, isoareag, fluxmaxg, ellipg
readcol,"/Users/jkrick/palomar/LFC/catalog/SExtractor.i.cat",o, xwcsi, ywcsi, xcenteri, ycenteri, fluxautoi, magautoi, magerri, fwhmi, isoareai, fluxmaxi, ellipi
readcol,"/Users/jkrick/palomar/LFC/catalog/SExtractor.u.cat",o, xwcsu, ywcsu, xcenteru, ycenteru, fluxautou, magautou, magerru, fwhmu, isoareau, fluxmaxu, ellipu

;test new zpts
magautou = magautou + 0.3
magautog = magautog + 0.1
magautor = magautor +0.35
magautoi = magautoi +0.35


numostars = 7417
sdssstar = replicate({star,umag:0D, gmag:0D,rmag:0D, imag:0D},numostars)
n = 0l
openr, lun, "/Users/jkrick/palomar/LFC/SDSS/sdss.large.csv", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF, lun, junk,run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z
    sdssstar[n] = {star, u,g, r, i}
    n = n + 1
ENDWHILE
close, lun
free_lun, lun

plot, sdssstar.gmag, sdssstar.gmag - sdssstar.rmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-r", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 0.55, thick = 3, color = colors.gray
plot, magautog, magautog - magautor,   thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-r", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
;oplot, findgen(25), findgen(25) - findgen(25) + 0.75, thick = 3, color = colors.gray
oplot, findgen(25), findgen(25) - findgen(25) + 0.55, thick = 3, color = colors.gray

plot, sdssstar.gmag, sdssstar.gmag - sdssstar.imag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-i", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 0.7, thick = 3, color = colors.gray
plot, magautog, magautog - magautoi,   thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-i", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
;oplot, findgen(25), findgen(25) - findgen(25) + 1.0, thick = 3, color = colors.gray
oplot, findgen(25), findgen(25) - findgen(25) + 0.7, thick = 3, color = colors.gray

plot, sdssstar.rmag, sdssstar.rmag - sdssstar.imag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="r magnitude", ytitle="r-i", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 0.25, thick = 3, color = colors.gray
plot, magautor, magautor - magautoi,   thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="r-i", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 0.25, thick = 3, color = colors.gray

plot, sdssstar.gmag, sdssstar.umag - sdssstar.gmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="u-g", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 1.25, thick = 3, color = colors.gray
plot, magautog, magautou - magautog,   thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="u-g", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
;oplot, findgen(25), findgen(25) - findgen(25) + 1.1, thick = 3, color = colors.gray
oplot, findgen(25), findgen(25) - findgen(25) + 1.25, thick = 3, color = colors.gray

plot, sdssstar.rmag, sdssstar.umag - sdssstar.rmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="r magnitude", ytitle="u-r", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 1.7, thick = 3, color = colors.gray
plot, magautor, magautou - magautor,   thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="r magnitude", ytitle="u-r", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
;oplot, findgen(25), findgen(25) - findgen(25) + 1.1, thick = 3, color = colors.gray
oplot, findgen(25), findgen(25) - findgen(25) + 1.7, thick = 3, color = colors.gray

plot, sdssstar.imag, sdssstar.umag - sdssstar.imag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="i magnitude", ytitle="u-i", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
oplot, findgen(25), findgen(25) - findgen(25) + 1.9, thick = 3, color = colors.gray
plot, magautoi, magautou - magautoi,   thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="i magnitude", ytitle="u-i", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
;oplot, findgen(25), findgen(25) - findgen(25) + 1.1, thick = 3, color = colors.gray
oplot, findgen(25), findgen(25) - findgen(25) + 1.9, thick = 3, color = colors.gray





ps_close, /noprint, /noid

undefine, sdssstar

end
