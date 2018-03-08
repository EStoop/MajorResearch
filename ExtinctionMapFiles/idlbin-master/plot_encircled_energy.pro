pro plot_encircled_energy, planetname
; this program plots the encircled energy as a function of aperture size for a staring mode data set
; the purpose of this code is to look for changes in the slope of that curve as a function of time


  print, systime()

  planetinfo = create_planetinfo()
  chname = planetinfo[planetname, 'chname']
  ra_ref = planetinfo[planetname, 'ra']
  dec_ref = planetinfo[planetname, 'dec']
  aorname = planetinfo[planetname, 'aorname']
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')


  restore, strcompress(dirname + planetname +'_ee_ch'+chname+'.sav')
  ap =  findgen(45) / 5 + 1.   ;this will make a smoother curve


;test output 
 ; print, 'bin_flux', bin_flux(0,*)


;plot fluxes in an annulus as a function of time
; I think I can just subtract two adjacent apertures to get an annulus.

;what about the central pixel
t2 = plot(findgen(n_elements(h))/850*4.0, bin_flux(*,0), '1r', xtitle = 'Time(hrs)', ytitle = 'Binned fluxes', name = '1 pixel aper')

annul3_1 = bin_flux(*,15) - bin_flux(*,7)
t = plot(findgen(n_elements(h))/850*4.0, annul3_1 +0.15, '1r',color = 'black', /overplot, name = '3pix - 1pix annulus')
l = legend(target = [t2, t], /data, /AUTO_TEXT_COLOR)

;Plot the binned curve of growths
  min = 500
  for p = min, 800 do begin
     if p eq min then begin
;        a = plot(ap, bin_flux(p,*), '1r', xtitle = 'Aperture size', ytitle = 'Encircled Energy') 
     endif
     if p gt 594 and p lt 670 then begin
;        a = plot(ap, bin_flux(p,*), '1b',/overplot)
     endif else begin
;        a = plot(ap, bin_flux(p,*), '1r', /overplot)
     endelse
     
  endfor


;-----------------------------
;generate binned plots from the stacked image
;------------------------------


restore,  '/Users/jkrick/irac_warm/pcrs_planets/hd7924/bin_stack.sav'

;the center pixel as a function of time
m = plot( findgen(850) / 850. * 4., bin_stack(15.0,15.0,*) / max(bin_stack(15.0,15.0,*)), '1r', xtitle = ' TIme(hrs)', ytitle = 'Normalized Flux')

;square annulus from 3 - 1 pixels
mean_an = fltarr(850)
for j = 0, 850 - 1 do begin
   a = total(bin_stack(13:14, 13:17,j) )
   b =total(bin_stack(15:17, 16:17,j))
   c =total(bin_stack(16:17, 13:15,j))
   d =total(bin_stack(15, 13:14,j))
   sumall = a + b + c + d
   mean_an(j) = sumall / 24.
   
endfor


;plotting as a function of time
m2 = plot(findgen(850)/850. *4., mean_an /max(mean_an) - 0.08, '1r',color = black, /overplot)

  print, systime()

end
