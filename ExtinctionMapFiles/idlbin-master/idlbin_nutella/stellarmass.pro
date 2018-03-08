pro stellarmass

restore, 'object.sav'
wgal = 10


z = object[wgal].zphot
;z = 3.17
print, "redshift", z

;put everything into microjanskies
uflux = 1E6*(10^((object[wgal].umaga - 8.926)/(-2.5)))
gflux = 1E6*(10^((object[wgal].gmaga - 8.926)/(-2.5)))
rflux = 1E6*(10^((object[wgal].rmaga - 8.926)/(-2.5)))
iflux = 1E6*(10^((object[wgal].imaga - 8.926)/(-2.5)))
unoise = 1E6*((10^((object[wgal].umaga - 8.926)/(-2.5))) - (10^((object[wgal].umaga +object[wgal].umagerra- 8.926)/(-2.5))))
gnoise = 1E6*((10^((object[wgal].gmaga - 8.926)/(-2.5))) - (10^((object[wgal].gmaga +object[wgal].gmagerra- 8.926)/(-2.5))))
rnoise = 1E6*((10^((object[wgal].rmaga - 8.926)/(-2.5))) - (10^((object[wgal].rmaga +object[wgal].rmagerra- 8.926)/(-2.5))))
inoise = 1E6*((10^((object[wgal].imaga - 8.926)/(-2.5))) - (10^((object[wgal].imaga +object[wgal].imagerra- 8.926)/(-2.5))))

obs = [uflux, gflux, rflux, iflux, object[wgal].irac1,  object[wgal].irac2,  object[wgal].irac3,  object[wgal].irac4] 
noise =[unoise,gnoise,rnoise,inoise,0.05*object[wgal].irac1,  0.05*object[wgal].irac2,  0.05*object[wgal].irac3,  0.05*object[wgal] .irac4] 

;what are the predicted values?
;read from bc03 output
;will want to do all this reading somewhere else, and make an array to then search

;x = [.3540,.4660,.6255,.7680,3.6000,4.5000,5.8,8.]               ;microns
x = [3540.,4660.,6255.,7680.,36000.,45000.,58000.,80000.]  ;angstroms

readcol,'/Users/jkrick/bin/bc03/models/Padova1994/salpeter/hr_m32_salp_ssp_age0.1,0.5,1,3,5,8,12',lambda, $
        exflux1,exflux5,exflux10,exflux30,exflux50,exflux80,exflux120,format="A"

;find the redshifted wavelengths in the model stellar population
point = fltarr(n_elements(x))
for i = 0, n_elements(x)-1 do begin
   if i lt 4 then pm = 9. else pm = 50.
   for l = 0, n_elements(lambda) - 1 do begin
  
       if lambda(l) lt x(i)/(1+z)+ pm and lambda(l) gt x(i)/(1+z) - pm then begin
         point(i) = l
         print, "found a match,x(i)lambda(l),l,i   ", x(i), " ",x(i)/(1+z), " ", lambda(l), exflux1(l)
      endif

   endfor
endfor

ex1 = exflux1(point)
ex5 =  exflux5(point)
ex10 = exflux10(point)
ex30 = exflux30(point)
ex50 = exflux50(point)
ex80 = exflux80(point)
ex120 = exflux120(point)


ex = [ex1,ex5,ex10,ex30,ex50,ex80,ex120]

;plot, lambda, ex1,/xlog
;oplot, lambda, exflux5
;oplot, lambda, exflux10
;oplot, lambda, exflux30
;oplot, lambda, exflux50
;oplot, lambda, exflux80
;oplot, lambda, exflux120

;convert from solar luminosities per angstrom(bc03 output) to microJy(sensible)
factor = 0.01266; 3.8E39 / 3E41;((3.826E33*1.E6)/(3.E18*1.E23))
ex = ex*factor

;ex = [1.09,1.42,6.57, 10.95, 34.52, 28.46, 26.72, 32.64]


;is there some way to scale the flux to make the chi-sqrd better?
start = [5E6]

result1= MPFITFUN('linear2',ex1,obs, noise, start,/quiet) 
result5= MPFITFUN('linear2',ex5,obs, noise, start,/quiet) 
result10= MPFITFUN('linear2',ex10,obs, noise, start,/quiet) 
result30= MPFITFUN('linear2',ex30,obs, noise, start,/quiet) 
result50= MPFITFUN('linear2',ex50,obs, noise, start,/quiet) 
result80= MPFITFUN('linear2',ex80,obs, noise, start,/quiet) 
result120= MPFITFUN('linear2',ex120,obs, noise, start,/quiet) 
;plot, obs, result(0)*obs + result(1), thick = 3,psym = 2;, color = colors.green


;scale the flux by the fitted value
;will need this value to scale the stellar mass
ex1 = result1(0)*ex1
ex5 = result5(0)*ex5
ex10 = result10(0)*ex10
ex30 = result30(0)*ex30
ex50 = result50(0)*ex50
ex80 = result80(0)*ex80
ex120 = result120(0)*ex120

;determine the chi squared of the fit between two arrays
df = n_elements(ex1) - 1  ;degrees of freedom

residual = (obs - ex120)
chi = total(residual^(2.0) / noise^(2.0))
prob = 1 - chisqr_pdf(chi, df)
print, chi, prob, result1(0)

plot, x, obs,thick =3,psym = 2,/xlog,/ylog, yrange=[0.1,100]
oplot, lambda*(1+z), exflux80*result80(0)
oplot, lambda*(1+z), exflux120*result120(0),linestyle = 2


end


;ex1 = ex1*factor
;ex5 = ex5*factor
;ex10 = ex10*factor
;ex30 = ex30*factor
;ex50 = ex50*factor
;ex80 = ex80*factor
;ex120 = ex120*factor


;ex1 = [exflux1(235),exflux1(1355),exflux1(2950),exflux1(4375),exflux1(6331),exflux1(6376),exflux1(6441),exflux1(6510)]
;ex5 = [exflux5(235),exflux5(1355),exflux5(2950),exflux5(4375),exflux5(6331),exflux5(6376),exflux5(6441),exflux5(6510)]
;ex10 = [exflux10(235),exflux10(1355),exflux10(2950),exflux10(4375),exflux10(6331),exflux10(6376),exflux10(6441),exflux10(6510)]
;ex30 = [exflux30(235),exflux30(1355),exflux30(2950),exflux30(4375),exflux30(6331),exflux30(6376),exflux30(6441),exflux30(6510)]
;ex50 = [exflux50(235),exflux50(1355),exflux50(2950),exflux50(4375),exflux50(6331),exflux50(6376),exflux50(6441),exflux50(6510)]
;ex80 = [exflux80(235),exflux80(1355),exflux80(2950),exflux80(4375),exflux80(6331),exflux80(6376),exflux80(6441),exflux80(6510)]
;ex120 = [exflux120(235),exflux120(1355),exflux120(2950),exflux120(4375),exflux120(6331),exflux120(6376),exflux120(6441),exflux120(6510)]
