pro pcrstiming

ps_start, filename= '/Users/jkrick/irac_warm/pcrsstare/timing.ps'
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

;start of the first pcrsstare
pcrsstare_start = 238.7778

;here are the results including x and y centroids from the pcrsstare observations
readcol, '/Users/jkrick/irac_warm/pcrsstare/v2_x_y_cen.txt', vw, day, reqkey, x, xerr, y, yerr, sclk0, sclk1
yesvw = where(vw eq 1)
novw = where(vw eq 0)
centroid_time = pcrsstare_start + day

for i = 0, n_elements(x) - 1 do print, centroid_time(i), x(i)

;here is the schedule summary of that week

readcol,'/Users/jkrick/irac_warm/pcrsstare/wk510_summary.csv', junk ,  aorkey , start,  title , target,format='(A,A,A,A,A)', delim = ',',/silent

;re-format start/stop keyword into floats 
start_days = float(strmid(start, 5, 3))
start_hrs = float(strmid(start, 9, 2))
start_mins = float(strmid(start, 12, 2))
start_secs = float(strmid(start, 15, 4))


;need to turn time into decimal time ie 12 hrs = 0.5 days
time = start_days + (start_hrs+( (start_mins + start_secs/60.)/60.)) / 24.

;make a plot
;st = plot(centroid_time(yesvw), x(yesvw),'g2D-',yrange = [14.3, 15.3], xrange = [244,253], xtitle = 'Time(days)',ytitle = 'X pix',  title = 'wk510 timeline')
;st = plot(centroid_time(novw), x(novw),'r2D-',/overplot)
 
 plot, centroid_time(yesvw), x(yesvw),yrange = [14.3, 15.3], xrange = [238,251], xtitle = 'Time(days)',ytitle = 'X pix',  title = 'wk510 timeline',  xstyle = 1,/nodata
oplot, centroid_time(yesvw), x(yesvw),color = greencolor

oplot,centroid_time(novw), x(novw), color = redcolor

line_x = [244,253]
line_y = [14.75, 14.75]
;line5 = polyline(line_x, line_y, thick = 1, color = red,/data)

line_x = [244,253]
line_y = [15.125, 15.125]
;line5 = polyline(line_x, line_y, thick = 1, color = green,/data)

;for a = 0,  n_elements(title) -1 do begin
;   xyouts, time(a),  14.3, title(a),orientation = 90, charsize = 1
;endfor

;new plan
b = where(title eq 'pcrsbias')
for count = 0,  n_elements(b) - 1 do xyouts, time(b), 14.3, title(b), color = bluecolor, orientation = 90

b = where(title eq 'pcrsstare')
for count = 0,  n_elements(b) - 1 do xyouts, time(b), 14.3, title(b), color = greencolor, orientation = 90

b = where(title eq 'pmap')
for count = 0,  n_elements(b) - 1 do xyouts, time(b), 14.3, title(b), color = purplecolor, orientation = 90

b = where(title eq 'iru_3axis_tng')
for count = 0,  n_elements(b) - 1 do xyouts, time(b), 14.3, title(b), color = redcolor, orientation = 90

b = where(title eq 'inc_c')
for count = 0,  n_elements(b) - 1 do xyouts, time(b), 14.3, title(b), color = cyancolor, orientation = 90

b = where(title eq 'pcrsbright' or title eq 'pcrsfaint' or title eq 'pcrsother')
for count = 0,  n_elements(b) - 1 do xyouts, time(b), 14.3, title(b), color = orangecolor, orientation = 90

al_legend, ['pcrsbias', 'pcrsstare', 'pmap', 'iru_3axis_tng', 'inc_c', 'bright/faint/other'], colors = ['blue', 'green', 'purple','red','cyan','orange'], position = [246, 15.3]

ps_end, /png

end
