function fisher_yates_shuffle, a
  n = n_elements(a) 
  if n lt 1 then begin
    print,' usage: fisher_yates_shuffle,SomeArray'
    print,'   Elements of SomeArray are placed in random order'
    return,0
  endif
  for i=n-1,1,-1 do begin
    ;; pick at random one of the first i elements
    indy = floor(randomu(seed)*(i+1))
    ;; move this element to the ith position
    tmp = a[i]
    a[i] = a[indy]
    a[indy] = tmp
 endfor
return, a
end
