pro read_ssc_query

readcol, '/Users/jkrick/external/irac_warm/darks/superdarks/iracmappc_subarray_2010thru2014.txt', pid, reqkey, duration, ch1, ch2, frametime, startday, starthour, aor_label, format = '(I10, I10, F10.2, A, A, F10.2, A, A, A)', comment = ';'
timestamptovalues, startday, year = year, month = month, day = day
xdata = [12, 4]
ydata = xdata

for y = 2010, 2013 do begin

;ch2 only, which also removes the darks
   t = where(year eq y and ch1 eq 'f' and ch2 eq 't')
   montht = month(t)

   plothist, montht, xhist, yhist, bin = 1, /noprint, /noplot
   print, 'yhist', yhist

   case y of 
      2010: begin
         ph10 = barplot(xhist, yhist, xtitle = 'Month', ytitle = 'Number of ch2 Sub Obs.', fill_color = [8,104,172], name = '2010')
         ydata2010 = yhist
         ;print, 'ydata2010', ydata2010
      end

      2011: begin
         yhist = [yhist, 0]
         xhist = findgen(12) + 1
         ydata2011 = ydata2010 + yhist
         ;print, 'ydata2011', ydata2011
         ph11 = barplot(xhist, ydata2011, bottom_values = ydata2010, fill_color = [67,162,202], /overplot, name = '2011')
      end
      2012:begin
         ydata2012 = ydata2011 + yhist
         print, 'ydata2012', ydata2012
         ph12 = barplot(xhist, ydata2012, bottom_values = ydata2011, fill_color = [123,204,196], /overplot, name = '2012')
      end
      2013:begin
         ydata2013 = ydata2012 + yhist
         print, 'ydata2012', ydata2013
         ph13 = barplot(xhist, ydata2013, bottom_values = ydata2012, fill_color = [186,228,188], /overplot, name = '2013')
      end

   endcase

endfor
   
leg = legend(target = [ph10,ph11, ph12, ph13], position = [4, 130], /data, /auto_text_color)
end
