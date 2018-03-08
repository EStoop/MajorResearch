function track_box, im, x_center, y_center, verbose = verbose
  pixvals = fltarr(49,63)
  for ip = 0, n_elements(x_center) - 1 do begin
     c = 0
     for xp = round(x_center(ip)) - 4, round(x_center(ip)) + 4 do begin
        for yp = round(y_center(ip)) - 4, round(y_center(ip)) + 4 do begin
           if keyword_set(verbose) then print, 'x,y', xp, yp,im[xp, yp]
           if xp gt 0 and yp gt 0 and xp lt 32 and yp lt 32 then pixvals[c] = im[xp, yp]
           c = c + 1
        endfor
     endfor

  endfor

return, pixvals
end
