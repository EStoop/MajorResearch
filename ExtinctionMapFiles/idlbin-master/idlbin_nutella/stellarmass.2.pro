pro stellarmass

;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/ZPHOT/stellarmass.2.ps", /portrait, xsize = 4, ysize = 4,/color

restore, '/Users/jkrick/idlbin/object.sav'


print,"started at "+systime()


;for this galaxy
for wgal = 0,4 do begin
   if object[wgal].prob ge 70 then begin
      ;don't bother calculating stellar mass without a good z calculation
   
      if object[wgal].gfwhm lt 5.6 and object[wgal].gfwhm gt 0 then begin
         ;don't calculate for stars
         print, "rejecting a star"
      endif else begin
      print, "galaxy ", wgal, object[wgal].zphot, object[wgal].jmag, object[wgal].tmassj, object[wgal].gfwhm

      z = object[wgal].zphot

       ;put everything into microjanskies
      ;don't include the non measurements, or bad measurements
      x = fltarr(12)   ;11 wavelengths
      obs = fltarr(12)
      noise = fltarr(12)
      
      nfluxes = 0
      if object[wgal].umaga lt 30 and object[wgal].umaga gt 10 then begin
         x(nfluxes) = 3540.
         obs(nfluxes) = 1E6*(10^((object[wgal].umaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((object[wgal].umaga - 8.926)/(-2.5))) - (10^((object[wgal].umaga +object[wgal].umagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes + 1
      endif
      if object[wgal].gmaga lt 30 and object[wgal].gmaga gt 10 then begin
         x(nfluxes) = 4660.
         obs(nfluxes) = 1E6*(10^((object[wgal].gmaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((object[wgal].gmaga - 8.926)/(-2.5))) - (10^((object[wgal].gmaga +object[wgal].gmagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes + 1
      endif
      if object[wgal].rmaga lt 30 and object[wgal].rmaga gt 10 then begin
         x(nfluxes) = 6255.
         obs(nfluxes) =1E6*(10^((object[wgal].rmaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((object[wgal].rmaga - 8.926)/(-2.5))) - (10^((object[wgal].rmaga +object[wgal].rmagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes + 1
      endif
      if object[wgal].imaga lt 30 and object[wgal].imaga gt 10 then begin
         x(nfluxes) = 7680.
         obs(nfluxes) =1E6*(10^((object[wgal].imaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((object[wgal].imaga - 8.926)/(-2.5))) - (10^((object[wgal].imaga +object[wgal].imagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes+1
      endif
      if object[wgal].jmag lt 30 and object[wgal].jmag gt 10 then begin
         x(nfluxes) = 12490.
         obs(nfluxes) =1E6*1594.*10^(object[wgal].jmag/(-2.5))
         noise(nfluxes) = 1E6*((10^((object[wgal].jmag - 8.926)/(-2.5))) - (10^((object[wgal].imaga +object[wgal].jmagerr- 8.926)/(-2.5))))
         nfluxes = nfluxes+1
      endif

     if object[wgal].tmassj ne 99 and object[wgal].tmassj gt 0 then begin
         x(nfluxes) = 12400.
         obs(nfluxes) = 1E6*object[wgal].tmassj
         noise(nfluxes) = 0.05*1E6*object[wgal].tmassj
         nfluxes = nfluxes+1
      endif
     if object[wgal].tmassh ne 99 and object[wgal].tmassh gt 0 then begin
         x(nfluxes) = 16500.
         obs(nfluxes) = 1E6*object[wgal].tmassh
         noise(nfluxes) = 0.05*1E6*object[wgal].tmassh
         nfluxes = nfluxes+1
      endif
     if object[wgal].tmassk ne 99 and object[wgal].tmassk gt 0 then begin
         x(nfluxes) = 21600.
         obs(nfluxes) = 1E6*object[wgal].tmassk
         noise(nfluxes) = 0.05*1E6*object[wgal].tmassk
         nfluxes = nfluxes+1
      endif
      if object[wgal].irac1 ne 99 and object[wgal].irac1 gt 0 then begin
         x(nfluxes) = 36000.
         obs(nfluxes) = object[wgal].irac1
         noise(nfluxes) = 0.05*object[wgal].irac1
         nfluxes = nfluxes+1
      endif
     if object[wgal].irac2 ne 99 and object[wgal].irac2 gt 0 then begin
         x(nfluxes) = 45000.
         obs(nfluxes) = object[wgal].irac2
         noise(nfluxes) = 0.05*object[wgal].irac2
         nfluxes = nfluxes+1
      endif
     if object[wgal].irac3 ne 99 and object[wgal].irac3 gt 0 then begin
         x(nfluxes) = 58000.
         obs(nfluxes) = object[wgal].irac3
         noise(nfluxes) = 0.05*object[wgal].irac3
         nfluxes = nfluxes+1
      endif
     if object[wgal].irac4 ne 99 and object[wgal].irac4 gt 0 then begin
         x(nfluxes) = 80000.
         obs(nfluxes) = object[wgal].irac4
         noise(nfluxes) = 0.05*object[wgal].irac4
         nfluxes = nfluxes+1
      endif

     x = x[0:nfluxes - 1]
     obs = obs[0:nfluxes-1]
     noise = noise[0:nfluxes-1]
     

;----------------------------------------------------------------------------
    ;what are the predicted values?

   ;read in the models
      readcol,'/Users/jkrick/bin/bc03/filelist',filename,format="A",/silent
      filename = strcompress(filename + ",0.5,1,3,5,8,12")

      finalchi = 1E6

   ;for each model
      for model = 0, n_elements(filename) -1 do begin

         readcol,filename(model),lambda,  exflux1,exflux5,exflux10,exflux30,exflux50,exflux80,exflux120,format="A",/silent
         exfluxar = [exflux1, exflux5,exflux10,exflux30,exflux50,exflux80,exflux120]

       ;find the redshifted wavelengths in the model stellar population
         point = fltarr(n_elements(x))
         for i = 0, n_elements(x)-1 do begin
            if i lt 4 then pm = 9. else pm = 50.
            for l = 0, n_elements(lambda) - 1 do begin
               if lambda(l) lt x(i)/(1+z)+ pm and lambda(l) gt x(i)/(1+z) - pm then begin
                  point(i) = l
               endif
            endfor
         endfor
 
      ;expected values as a function of age, at the correct wavelengths
      ;for comparison with the observed values
         ex1 = exflux1(point)
         ex5 =  exflux5(point)
         ex10 = exflux10(point)
         ex30 = exflux30(point)
         ex50 = exflux50(point)
         ex80 = exflux80(point)
         ex120 = exflux120(point)

        
         ex = [ex1,ex5,ex10,ex30,ex50,ex80,ex120]
      
  
      ;convert from solar luminosities per angstrom(bc03 output) to microJy(sensible)
         factor = 0.01266       ; 3.8E39 / 3E41;((3.826E33*1.E6)/(3.E18*1.E23))
         ex = ex*factor

         ex1 = ex1*factor
         ex5 = ex5*factor
         ex10 = ex10*factor
         ex30 = ex30*factor
         ex50 = ex50*factor
         ex80 = ex80*factor
         ex120 = ex120*factor

      ;scale the flux to make the chi-sqrd better?
         start = [5E6]

         result1= MPFITFUN('linear2',ex1,obs, noise, start,/quiet) 
         result5= MPFITFUN('linear2',ex5,obs, noise, start,/quiet) 
         result10= MPFITFUN('linear2',ex10,obs, noise, start,/quiet) 
         result30= MPFITFUN('linear2',ex30,obs, noise, start,/quiet) 
         result50= MPFITFUN('linear2',ex50,obs, noise, start,/quiet) 
         result80= MPFITFUN('linear2',ex80,obs, noise, start,/quiet) 
         result120= MPFITFUN('linear2',ex120,obs, noise, start,/quiet) 
      ;plot, ob s, result(0)*obs + result(1), thick = 3,psym = 2;, color = colors.green

         result = [result1, result5,result10,result30,result50,result80,result120]

      ;scale the flux by the fitted value
      ;will need this value to scale the stellar mass
         ex1 = result1(0)*ex1
         ex5 = result5(0)*ex5
         ex10 = result10(0)*ex10
         ex30 = result30(0)*ex30
         ex50 = result50(0)*ex50
         ex80 = result80(0)*ex80
         ex120 = result120(0)*ex120

         extest = obs*0.9
         ex = [ex1,ex5,ex10,ex30,ex50,ex80,ex120]




      ;determine the chi squared of the fit between two arrays
         df = n_elements(ex1) - 1 ;degrees of freedom

         for nage = 0, n_elements(result) - 1 do begin
            b = nage*n_elements(x)
            e = nage*n_elements(x)+ n_elements(x) - 1
            residual = (obs - ex(b:e))
            chi = total(residual^(2.0) / noise^(2.0))
            prob = 1 - chisqr_pdf(chi, df)
;            print, "age, chi, prob, result(0)", nage, chi, prob ;, result(nage)

            if chi lt finalchi then begin
               finalchi = chi
               finalprob = prob
               mass = result(nage)
              print, "best chi so far", chi, prob
               b = nage*1085
               e = nage*1085+1084
              plot, x, obs, psym = 2, thick = 3,title = wgal, color = colors.red,/xlog,/ylog
               oplot, lambda*(1+z), exfluxar(b:e)*factor*result(nage)
 ;              oplot, lambda, exfluxar(b:e)*factor*result(nage), linestyle =2
               
            endif

         endfor



      endfor
      print, "smallest chi", finalchi, finalprob, mass
      endelse
   endif
endfor


print,"Finished at "+systime()

ps_close, /noprint, /noid
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

;plot, lambda, ex1,/xlog
;oplot, lambda, exflux5
;oplot, lambda, exflux10
;oplot, lambda, exflux30
;oplot, lambda, exflux50
;oplot, lambda, exflux80
;oplot, lambda, exflux120

;        oplot, lambda*(1+z), exflux5*factor*result5(0),linestyle = 1
;         oplot, lambda*(1+z), exflux10*factor*result10(0),linestyle = 2
;         oplot, lambda*(1+z), exflux30*factor*result30(0),linestyle = 3
;         oplot, lambda*(1+z), exflux50*factor*result50(0),linestyle = 4
;         oplot, lambda*(1+z), exflux80*factor*result80(0),linestyle = 5
;         oplot, lambda*(1+z), exflux120*factor*result120(0)

;names = ["exflux1", "exflux5", "exflux10","exflux30","exflux50","exflux80","exflux120"]


;      obs = [uflux, gflux, rflux, iflux, object[wgal].irac1,  object[wgal].irac2,  object[wgal].irac3,  object[wgal].irac4] 
;     noise =[unoise,gnoise,rnoise,inoise,0.05*object[wgal].irac1,  0.05*object[wgal].irac2,  0.05*object[wgal].irac3,  0.05*object[wgal] .irac4] 
   ;x = [.3540,.4660,.6255,.7680,3.6000,4.5000,5.8,8.]               ;microns
;      x = [3540.,4660.,6255.,7680.,36000.,45000.,58000.,80000.] ;angstroms
