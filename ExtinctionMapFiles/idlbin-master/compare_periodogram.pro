pro compare_periodogram

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/periodogram1.sav'
;b2 = plot(wk1, wk2, /overplot, color = 'black',thick = 2, name = 'Original')

b = plot(wk1, wk2, xtitle = 'Frequency(1/days)', ytitle = 'Power', xrange =[0,50], yrange = [0,150],thick = 2, color = 'aqua',name = 'Pmap Corrected');

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/periodogram2.sav'
b2 = plot(wk1, wk2, /overplot, color = 'blue',thick = 2, name = 'blue')

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/periodogram3.sav'
b2 = plot(wk1, wk2, /overplot, color = 'steel_blue',thick = 2, name = 'blue')

b2.save, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/compare_periodograms.png'
end
