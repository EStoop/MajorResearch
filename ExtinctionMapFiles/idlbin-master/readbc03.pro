pro readbc03

redcolor = FSC_COLOR("Red", !D.Table_Size-2);
orangecolor = FSC_COLOR("Orange", !D.Table_Size-3);
greencolor = FSC_COLOR("Green", !D.Table_Size-4);
bluecolor = FSC_COLOR("Blue", !D.Table_Size-5)
purplecolor = FSC_COLOR("Purple", !D.Table_Size-6)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-7)

filename = '/Users/jkrick/bc03/models/Padova1994/chabrier/bc2003_lr_m32_chab_ssp.spect'

readcol, filename, lambda, age0, age1, age2,age3,age4,age5,age6,format="A",/silent

plot, lambda, age0/ age0(989), charthick = 1, /ylog,/xlog, thick = 1, xrange = [3E4, 6E4], xstyle = 1, yrange = [.1,10.]
oplot, lambda, age1/ age1(989), thick = 1, color = purplecolor
oplot, lambda, age2/ age2(989), thick = 1, color = bluecolor
oplot, lambda, age3/ age3(989), thick = 1, color = greencolor
oplot, lambda, age4/ age4(989), thick = 1, color = yellowcolor
oplot, lambda, age5/ age5(989), thick = 1, color = orangecolor
oplot, lambda, age6/ age6(989), thick = 1, color = redcolor



end
