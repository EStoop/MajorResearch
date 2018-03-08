pro test_nearest_neighbor_array
x = randomu(1, 100)*100
y = randomu(2, 100)*100
sqrtnp = randomu(3, 100)*100

dn = nth_neighbor(x, y, 2)
help, dn
for i = 0, n_elements(x) - 1 do begin
   print, x(i), y(i), x(dn(i)), y(dn(i))
endfor

end

function nth_neighbor, x, y, k
  n = n_elements(x)
  dn = fltarr(n,/NOZERO)
  d=(rebin(transpose(x),n,n,/SAMPLE)-rebin(x,n,n,/SAMPLE))^2 + $
    (rebin(transpose(y),n,n,/SAMPLE)-rebin(y,n,n,/SAMPLE))^2 
  for i=0L,n-1 do dn[i] = sqrt(d[(sort(d[*,i]))[k],i])
  return, dn
end
