pro mass_test
  close,/all
  !p.multi = [0, 0, 1]
  ps_open, file = "/Users/jkrick/maraston/compare_mass.ps", /portrait, xsize = 4, ysize = 4,/color
 
 restore, '/Users/jkrick/idlbin/object.maraston.sav'
  marmass = object.mass
  
  print, marmass(0:30)
  restore, '/Users/jkrick/idlbin/object.sav'
  bcmass = object.mass
  print, bcmass(0:30)
  plot, alog10(marmass), alog10(bcmass),psym = 2, xrange=[4,14],yrange=[4,14],xstyle = 1, ystyle = 1,$
        xtitle = "maraston mass", ytitle = "bc03 mass"
  
  oplot, findgen(20), findgen(20)
  ps_close, /noprint, /noid
end
