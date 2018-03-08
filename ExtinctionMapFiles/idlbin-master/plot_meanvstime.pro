pro plot_meanvstime
restore,  '/Users/jkrick/virgo/irac/fitmean_ch1.sav'

ps_start, filename= '/Users/jkrick/virgo/irac/meanvstime_ch1.ps',/nomatch

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
yellowcolor = FSC_COLOR("yellow", !D.Table_Size-9)
colorarr = [redcolor, orangecolor, yellowcolor, greencolor, cyancolor, bluecolor, purplecolor]

!P.multi = [0,1,1]
!P.Thick = 3
!P.CharThick = 3
!P.Charsize = 1.25 
!X.Thick = 3
!Y.Thick = 3

plot, xt, meanarr,  xtitle = 'Frame Number', ytitle = 'Mean Levels' ,yrange = [0.22, 0.28]; yrange = [0.0,0.2];, color = colorarr[a] 

ps_end
end
