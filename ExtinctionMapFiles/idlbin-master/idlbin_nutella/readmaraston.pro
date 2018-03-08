pro readmaraston
close, /all
nlambda = 1220.D
;nages, nlambdas, nmodels
model = fltarr(67,1221,68)

readcol,'/Users/jkrick/maraston/Sed_Mar05_SSP_Salpeter/filelistlong', filename,format="A",/silent


for i = 0, n_elements(filename) -1 do begin

   readcol, filename(i), age, metal,  lambda, flux,format="A";,/silent
;   readcol, filename(i), lambda, age0, age1, age2,age3,age4,age5,age6,age7,age8,format="A",/silent
   
   print, "i, filename", i, filename(i)
   print, "n_elements(flux)", n_elements(flux)

   for j = 0, uint(n_elements(flux) / 1220) - 1 do begin
      model[j,*,i] = flux(nlambda*j +j : nlambda*(j+1) + j)

   endfor

endfor
;------------------------------------------------------------------
readcol,'/Users/jkrick/maraston/Sed_Mar05_SSP_Salpeter/filelistshort', filename,format="A",/silent


for i = i, n_elements(filename) +i -1 do begin
   n = 0
   readcol, filename(n), age, metal,  lambda, flux,format="A";,/silent
;   readcol, filename(i), lambda, age0, age1, age2,age3,age4,age5,age6,age7,age8,format="A",/silent
   
   print, "i, filename", i, filename(n)
   print, "n_elements(flux)", n_elements(flux)

   for j = 51, 67 - 1 do begin
 
      model[j,*,i] = flux(nlambda*n +n : nlambda*(n+1) + n)
     n = n + 1
 
   endfor
endfor
;------------------------------------------------------------------

readcol,'/Users/jkrick/maraston/Sed_Mar05_SSP_Salpeter/filelistveryshort', filename,format="A",/silent
;10 and  15 gyr  ie 61 and 66

for i = i, n_elements(filename)+i -1 do begin
   n = 0
   readcol, filename(n), age, metal,  lambda, flux,format="A";,/silent
;   readcol, filename(i), lambda, age0, age1, age2,age3,age4,age5,age6,age7,age8,format="A",/silent
   
   print, "i, filename", i, filename(n)
   print, "n_elements(flux)", n_elements(flux)

   for j = 61, 66,5 do begin
      model[j,*,i] = flux(nlambda*n +n : nlambda*(n+1) + n)
      n=n+1
   endfor

endfor

;csp----------------------------------------------------------------
readcol,'/Users/jkrick/maraston/filelist', filename,format="A",/silent

for i = i, n_elements(filename) +i -1 do begin
   n = 0
   readcol, filename(n), age, lambda, flux,format="A";,/silent
;   readcol, filename(i), lambda, age0, age1, age2,age3,age4,age5,age6,age7,age8,format="A",/silent
   
   print, "i, filename", i, filename(n)
   print, "n_elements(flux)", n_elements(flux)

   for j = 6, 67 - 1 do begin
      model[j,*,i] = flux(nlambda*n +n : nlambda*(n+1) + n)
      n = n + 1
   endfor

endfor
;------------------------------------------------------------------
lambda = lambda[0:1220]    ;just one pass through.
mmodel = model
save, mmodel, filename='/Users/jkrick/maraston/model.2.sav'
mlambda =lambda
save, mlambda, filename = '/Users/jkrick/maraston/lambda.sav'

readcol,'/Users/jkrick/maraston/agegridSSP_Mar05.txt', age,format="A",/silent
;save, age, filename = '/Users/jkrick/maraston/age.sav'

end
