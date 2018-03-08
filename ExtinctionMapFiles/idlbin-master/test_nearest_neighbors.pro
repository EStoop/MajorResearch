pro test_nearest_neighbors
x = randomu(sd, 100)*100
y = randomu(sd, 100)*100
sqrtnp = randomu(sd, 100)*100

nearest = nearest_neighbors_DT(x,y,'2', DISTANCES=nearest_d,NUMBER=3)
for i = 0, n_elements(x) - 1 do begin
   print, i, x(i), y(i), x(nearest(0,i)), y(nearest(0,i)),x(nearest(1,i)), y(nearest(1,i)), x(nearest(2,i)), y(nearest(2,i))
endfor

print, '-------------------------------'

nearest_np = nearest_neighbors_np_DT(x,y,sqrtnp,'2', DISTANCES=nearest_d,NUMBER=3)
for i = 0, n_elements(x) - 1 do begin
   print, i, x(i), y(i), x(nearest(0,i)), y(nearest(0,i)),x(nearest(1,i)), y(nearest(1,i)), x(nearest(2,i)), y(nearest(2,i))
endfor
end
