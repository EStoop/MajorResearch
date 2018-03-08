pro test_match2D
x1 = findgen(10)
y1 = findgen(10)
;print, x1
x2 = randomu(sd, 100)*100
y2 = randomu(sd, 100)*100
;print, x2
 match = match_2D(x1, y1, x2, y2, 50., match_distance = match_d)
 help, match
for i = 0, n_elements(x1) - 1 do begin
   if match(i) gt 0 then print, 'match, i ',i, x2[match(i)], y2[match(i)]
endfor

end
