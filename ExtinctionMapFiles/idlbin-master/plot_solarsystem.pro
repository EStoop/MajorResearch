pro plot_solarsystem

  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/comet-orbital-trends-colors.txt', name, rh, color, colorerr, class, format = '(A, D10, D10, D10, A)'
  readcol,  '/Users/jkrick/external/irac_warm/senior_review_2016/comet-orbital-trends-dust-templates.txt', rh_model, two, three, four, format = '(D10, D10, D10, D10)'
  SP = where(class eq 'SP', complement = OC)

  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/reach_comets.txt', name,Date,rh_reach,	dS,Solar_elong,	V_a, Fch1,Fch2,	Apertureb, format = '(A, A, F10, F10, F10, F10,F10,F10,F10)'
  
  p1 = errorplot(rh(SP), color(SP), colorerr(SP), '1s', xtitle = 'Distance from the Sun (AU)', ytitle = '[3.6] - [4.5] (mag)', sym_filled = 1, color = 'blue', yrange = [-1, 4], xrange = [0, 5], name = 'Short Period', errorbar_color = 'blue', xthick = 2, ythick = 2, xtickfont_style = 1, ytickfont_style = 1)
  p3 = errorplot(rh(OC), color(OC), colorerr(OC), '1o', color = 'red', sym_filled = 1, overplot = p1,errorbar_color =  'red', name = 'Oort Cloud')
  p2 = plot(rh_model, two, color = 'black', thick = 2, overplot = p1)
  p2 = plot(rh_model, three, color = 'black', thick = 2, overplot = p1)
  p2 = plot(rh_model, four, color = 'black', thick = 2, overplot = p1, name= 'Dust only models')
;  p4 = plot(rh_reach, Fch2 / Fch1, '1o', sym_size = 0.5, sym_filled = 1, color = 'black', overplot = p1)

;  t1 = text(3, 0, 'dust-only models', color = 'green', /data, overplot = p1)

  l = legend(target = [p1, p3, p2], position = [2.0, 3.8], /data, /auto_text_color)
end
