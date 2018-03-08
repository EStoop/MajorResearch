pro pcrsstare_v2

readcol, '/Users/jkrick/irac_warm/pcrsstare/v2_x_y_cen.txt', vw, day, reqkey, x, xerr, y, yerr, sclk0, sclk1
deltax = fltarr(n_elements(vw)/2.)
deltay = fltarr(n_elements(vw)/2.)
deltatime = fltarr(n_elements(vw)/2.)
a = 0
for i = 0, n_elements(vw) - 2, 2 do begin
   deltax[a] = x(i) - x(i + 1)
   deltay[a] = y(i) - y(i + 1)
   deltatime[a] = day(i)
   a = a + 1
endfor

s = plot(deltatime, abs(deltax), 'r2D-', xtitle = 'Days', ytitle = 'Delta x')
t = plot(deltatime, abs(deltay), 'b2D-', xtitle = 'Days', ytitle = 'Delta y')
end
