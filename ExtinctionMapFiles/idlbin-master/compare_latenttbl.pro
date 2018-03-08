pro compare_latenttbl

generatedname = '/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/latent_generated.tbl'
readcol, generatedname, aorname, dceid ,x, y,  currentsclk ,  triggersclk ,  size, format = '(A, L10, F10.2, F10.2,D26.8,D26.8, I10)'

;sort so that all flags on each dceid are listed together
a = sort(dceid)
dceid = dceid[a]
aorname = aorname[a]
x = x[a]
y =y[a]
currentsclk = currentsclk[a]
triggersclk  = triggersclk[a]

;openw, outlun, '/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/latent_generated_sort.tbl',/get_lun
for i = 0, n_elements(dceid) - 1 do begin ;
;   printf, outlun, strmid(aorname[i], 1,8), dceid[i], x[i], y[i], currentsclk[i], triggersclk[i], 10, format = '(A, I10, F10.2, F10.2,D26.8,D26.8, I10)'
endfor

;close, outlun
;free_lun, outlun

sortx = x(sort(x))
sorty = y(sort(x))

uniqx = sortx[uniq(sortx)]
uniqy = sorty[uniq(sortx)]

print, 'uniqx', n_elements(uniqx), uniqx

;-------------------------------

pipelinetblname =  '/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/latentflag.tbl'
readcol, pipelinetblname, paorname, pdceid ,px, py,  pcurrentsclk ,  ptriggersclk ,  psize, format = '(I10, L10, F10.2, F10.2,D26.8,D26.8, I10)'

sortpx = px(sort(px))
sortpy = py(sort(px))

uniqpx = sortpx[uniq(sortpx)]
uniqpy = sortpy[uniq(sortpx)]

print, 'uniqpx', n_elements(uniqpx), uniqpx


;-------------------
;make sure that every flagged pixel in my generated latent table is in the pipeline table
;and switch...

for i = 0, n_elements(uniqx) - 1 do begin
   a = where(uniqx(i) eq uniqpx, count)
;   if count lt 1 then print, 'no match', uniqx(i), uniqy(i)
endfor

for i = 0, n_elements(uniqpx) - 1 do begin
   a = where(uniqpx(i) eq uniqx, count)
   if count lt 1 then print, 'no match', uniqpx(i), uniqpy(i)
endfor

;---------------------------------------
;check that filling in the gaps happened.  Pick one trigger, and print out each dceid where that trigger occurs. 
;Are there any gaps (not HDR induced).

xcen =  224.26  
ycen =  222.820

a = where(xcen eq px and ycen eq py)
print, 'matched', paorname(a(1)),pdceid(a)

end
