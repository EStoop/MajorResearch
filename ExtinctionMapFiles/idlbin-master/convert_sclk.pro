pro convert_sclk

;want this number in sclk
;2011-247/23:09:31.4


;sclk is defined as the number of seconds since 1/1/1980
year = 2011
day = 247
hour = 23
min = 09
sec = 31.4

sclk = double(0)
sclk = sec + (min*60.) + (hour*60.*60.) + ( day*24.*60.*60.) + ((year - 1980)*365.*24.*60.*60.)

print, sclk, format = '(D14.2)'


; now go the other way
print, 'outliers'
outlier_sclk = [999405565.D,999781845.D,999904668.D,1000587195.D,1000590410.D]
outlier_sclk_ch2 = [998852449.D   ,998905550.D   ,999018147.D   ,999085582.D   ,999085967.D   ,999134696.D   ,999135082.D   ,999139537.D   ,999139922.D   ,999141794.D   ,999142180.D   ,999183786.D   ,999184171.D   ,999185144.D   ,999185531.D   ,999265374.D   ,999400218.D   ,999444386.D   ,999529793.D   ,999444987.D   ,999541376.D   ,999552600.D   ,999651532.D   , 999696224.D  , 1000552456.D ,1000589483.D  ,1000489782.D  ,1000588128.D  ]

outlier_sclk_2 = 939437952.D
outlier_sclk = outlier_sclk_ch2
dayarr = fltarr(n_elements(outlier_sclk))
for i = 0, n_elements(outlier_sclk) - 1 do begin
   year = (outlier_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   dayarr[i] = day
   print, sclk, year, day
endfor
youtlier = fltarr(n_elements(outlier_sclk)) + 1.
a = plot(dayarr, youtlier, '*r2', xtitle = 'day of year')

print, '-----------------------'
print, 'good'
good_sclk = [ 999444230.D ,  999693661.D ,  999732318.D ,  999751864.D ,  999756044.D ,  999760507.D ,  999769922.D ,  999771702.D ,  999546815.D ,  999553526.D ,  999577309.D ,  999644194.D ,  999686640.D ,  999689118.D ,  999782053.D , 1000072460.D ,  1000440302.D ,  999380117.D ,  999630508.D ,  999685707.D ,  999755111.D ,  999759574.D ,  999770131.D ,  1000257191.D , 1000205883.D ]
good_sclk_ch2 = [998852063.D   ,998905935.D   ,999018534.D   ,999217045.D   ,999217431.D   ,999265760.D   ,999272652.D   ,999273037.D   ,999382085.D   ,999382471.D   ,999527888.D   ,999528273.D   ,   999379515.D   ,   999391143.D   ,   999395786.D   ,   999405722.D   ,   999452723.D   ,   999540774.D   ,   999542758.D   ,   999688191.D   ,  1000204956.D   ,  1000258124.D   ,   999543360.D   ,   999545886.D   ,   999582183.D   ,   999802515.D   ,   999905598.D   , 1000159590.D   , 1000207404.D   ,  1000439372.D   ]

good_sclk = good_sclk_ch2
dayarr = fltarr(n_elements(good_sclk))

for i = 0, n_elements(good_sclk) - 1 do begin
   year = (good_sclk[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   dayarr[i]  = day
   print, year, day
endfor
ygood = fltarr(n_elements(good_sclk)) + 1.

a = plot(dayarr, ygood, '*g2',/overplot)

end



