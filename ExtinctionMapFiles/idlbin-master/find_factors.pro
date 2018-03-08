function find_factors, nsub, goodval
  
                                ; find all the factors of a number, and return the one closest to some preferred value

  factor, nsub, p, n
  factorarr = lonarr(fix(sqrt(nsub) / 2))
  ;want to find the one that is closest to 100
  c = 0
  for a = 0, n_elements(p) - 2 do begin
     for b = 0, n(a) do begin
        ;print,  (p(a))^b
        factorarr(c) =  (p(a))^b
        c = c + 1
        ;print,  p(a)^b * p(a+1)
        factorarr(c) =  p(a)^b * p(a+1)
        c = c + 1
     endfor
  endfor

  factorarr = factorarr[0:c-1]
  u = uniq(factorarr)
  factorarr = factorarr(u)
  s = sort(factorarr)
  factorarr = factorarr(s)
  ;print, factorarr
  

; find the one closest to goodval
; goodval = 60
test = abs(factorarr - goodval)
theone = where(test eq min(test))

return,  factorarr(theone)

end

