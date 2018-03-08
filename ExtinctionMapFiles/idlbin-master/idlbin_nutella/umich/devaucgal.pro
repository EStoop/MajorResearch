PRO devaucgal

close,/all

P0 = 21.5
P1 = 10
x = findgen(100)
yarr = fltarr(1000)
yarr =  P0 * (exp(-7.67*(((x/P1)^(1.0/4.0)) - 1.0)))

plot, x, yarr, thick = 3
END
