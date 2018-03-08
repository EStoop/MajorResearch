pro plot_superdark_seasonal
  dirname = '/Users/jkrick/external/irac_warm/darks/superdarks/'
  fits_read,  strcompress(dirname+ 'superdark_ch2_s2s.fits',/remove_all), superdark, header 

  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_feb.fits',/remove_all), superdark_feb, header 
  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_may.fits',/remove_all), superdark_may, header 
  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_aug.fits',/remove_all), superdark_aug, header ;
  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_nov.fits',/remove_all), superdark_nov, header 


;  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_2010.fits',/remove_all), superdark_feb, header 
;  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_2011.fits',/remove_all), superdark_may, header 
;  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_2012.fits',/remove_all), superdark_aug, header 
;  fits_read, strcompress(dirname+ 'superdark_ch2_s2s_2013.fits',/remove_all), superdark_nov, header 

  
;;--------------------------------------------------------------------------------
;what if I median together all 64 subframes
  superdark = median(superdark, dimension = 3)
  superdark_feb = median(superdark_feb, dimension = 3)
  superdark_may = median(superdark_may, dimension = 3)
  superdark_aug = median(superdark_aug, dimension = 3)
  superdark_nov = median(superdark_nov, dimension = 3)
;  superdark =superdark[*,*,0]
;  superdark_feb = superdark_feb[*,*,0]
;  superdark_may = superdark_may[*,*,0]
;  superdark_aug = superdark_aug[*,*,0]
;  superdark_nov = superdark_nov[*,*,0]

  superdark_feb = superdark - superdark_feb
  superdark_may = superdark - superdark_may
  superdark_aug = superdark - superdark_aug
  superdark_nov = superdark - superdark_nov


;;--------------------------------------------------------------------------------

  ;;try aperture photometry at the center and a lots of other random
  ;;locations
  
  ;;first find the random centers
  rx = (21*randomu(seed, 100)) + 5.
  ry =(21*randomu(seed, 100)) + 5.
  colorarr = ['sky_blue', 'blue', 'red',  'yellow']
  for n = 0, 3 do begin
     case n of
        0: data = superdark - superdark_feb
        1: data = superdark - superdark_may
        2: data = superdark - superdark_aug
        3: data = superdark - superdark_nov

      endcase

     aper, data, rx, ry, flux, fluxerr, bkgd, bkgderr, 1.0, 2.25, [3,7], $
           /FLUX, /EXACT, /NAN, /SILENT, /MEANBACK
     
     plothist, flux, xhist, yhist, bin = 1.0, /noplot, /noprint
     ph = barplot(xhist, yhist,  xtitle = 'Aperture Photometry', ytitle = 'Number', fill_color = colorarr[n], $
                  nbars = 4, index = n, overplot = ph, xrange = [-10, 10], yrange = [0,40], title = 'Yearly' )
     
     aper, data, 15.0, 15.0, flux, fluxerr, bkgd, bkgderr, 1.0, 2.25, [3,7], $
           /FLUX, /EXACT, /NAN, /SILENT, /MEANBACK
     print, 'flux at center', flux
     ph = plot([flux, flux], [0,40], linestyle = 2, thick = 2, overplot = ph, color = colorarr[n])
  endfor


;;print out the 64 stack
  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_feb_x64.fits',/remove_all), superdark_feb, header 
  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_may_x64.fits',/remove_all), superdark_may, header 
  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_aug_x64.fits',/remove_all), superdark_aug, header 
  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_nov_x64.fits',/remove_all), superdark_nov, header 
;  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_2010_x64.fits',/remove_all), superdark_feb, header 
;  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_2011_x64.fits',/remove_all), superdark_may, header 
;  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_2012_x64.fits',/remove_all), superdark_aug, header 
;  fits_write, strcompress(dirname+ 'diffdark_ch2_s2s_2013_x64.fits',/remove_all), superdark_nov, header 

end


;;quantile quantile plots 
 ;simplify for now
;  data1 = superdark_feb[*,*,0]
;  data2 = superdark_nov[*,*,0]

;  quant1 = quantile(data1, 20)
;  quant2 = quantile(data2, 20)

;  pq = plot(quant1, quant2, '1o', sym_filled = 1, xrange = [25, 50], yrange = [25, 50], color = 'red')
;  pq = plot(findgen(50), findgen(50), overplot = pq)
