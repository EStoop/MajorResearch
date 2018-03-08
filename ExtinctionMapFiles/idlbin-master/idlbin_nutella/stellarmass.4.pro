pro stellarmass

device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/ZPHOT/stellarmass.noz.ps", /portrait, xsize = 4, ysize = 4,/color

restore, '/Users/jkrick/idlbin/object.sav'
restore,'/Users/jkrick/bin/bc03/bc03.sav'
restore,'/Users/jkrick/bin/bc03/lambda.sav'

print,"started at "+systime()

readcol,'/Users/jkrick/bin/bc03/filelist',filename,format="A",/silent
filename = strcompress(filename + ",0.3,0.5,0.7,1,3,5,8,12")

;age = fltarr(7)
age = [0.1,0.3,0.5,0.7,1.,3.,5.,8.,12.]
;for this galaxy
for wgal = 0, 50 do begin; n_elements(object.rmaga) - 1 do begin
   if object[wgal].prob ge 70 then begin
      ;don't bother calculating stellar mass without a good z calculation
   
      if object[wgal].gfwhm lt 5.6 and object[wgal].gfwhm gt 0 then begin
         ;don't calculate for stars
         print, "rejecting a star"
      endif else begin
      print, "galaxy ", wgal, object[wgal].zphot

      z = object[wgal].zphot
      maxage = galage(z,1000,/silent)
      lumindist = lumdist(z,/silent)   ;in Mpc
      lumindist = lumindist * 3.08E24  ;in cm
      area = 0.D
      area = 4*!PI*(lumindist)^2     ;in cm^2
 
       ;put everything into microjanskies
      ;don't include the non measurements, or bad measurements
      x = fltarr(15)   ;15 wavelengths
      obs = fltarr(15)
      noise = fltarr(15)
      
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
      if object[wgal].flamjmag lt 30 and object[wgal].flamjmag gt 10 then begin
         x(nfluxes) = 12490.
         obs(nfluxes) =1E6*1594.*10^(object[wgal].flamjmag/(-2.5))
         noise(nfluxes) = 1E6*1594.*(10^(object[wgal].flamjmag/(-2.5)) - (10^((object[wgal].flamjmag+ object[wgal].flamjmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif
      if object[wgal].wircjmag lt 30 and object[wgal].wircjmag gt 10 then begin
         x(nfluxes) = 12500
         obs(nfluxes) =1E6*1594.*10^(object[wgal].wircjmag/(-2.5))
         noise(nfluxes) = 1E6*1594.*(10^(object[wgal].wircjmag/(-2.5)) - (10^((object[wgal].wircjmag+ object[wgal].wircjmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif

     if object[wgal].tmassj ne 99 and object[wgal].tmassj gt 0 then begin
         x(nfluxes) = 12400.
         obs(nfluxes) = 1E6*object[wgal].tmassj
         noise(nfluxes) = 0.05*1E6*object[wgal].tmassj
         nfluxes = nfluxes+1
      endif
      if object[wgal].wirchmag lt 30 and object[wgal].wirchmag gt 10 then begin
         x(nfluxes) = 16350
         obs(nfluxes) =1E6*1024.*10^(object[wgal].wirchmag/(-2.5))
         noise(nfluxes) = 1E6*1024.*(10^(object[wgal].wirchmag/(-2.5)) - (10^((object[wgal].wirchmag+ object[wgal].wirchmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif
     if object[wgal].tmassh ne 99 and object[wgal].tmassh gt 0 then begin
         x(nfluxes) = 16500.
         obs(nfluxes) = 1E6*object[wgal].tmassh
         noise(nfluxes) = 0.05*1E6*object[wgal].tmassh
         nfluxes = nfluxes+1
      endif
      if object[wgal].wirckmag lt 30 and object[wgal].wirckmag gt 10 then begin
         x(nfluxes) = 21500
         obs(nfluxes) =1E6*666.8*10^(object[wgal].wirckmag/(-2.5))
         noise(nfluxes) = 1E6*666.8*(10^(object[wgal].wirckmag/(-2.5)) - (10^((object[wgal].wirckmag+ object[wgal].wirckmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif
     if object[wgal].tmassk ne 99 and object[wgal].tmassk gt 0 then begin
         x(nfluxes) = 21600.
         obs(nfluxes) = 1E6*object[wgal].tmassk
         noise(nfluxes) = 0.05*1E6*object[wgal].tmassk
         nfluxes = nfluxes+1
      endif
      if object[wgal].irac1 ne 99 and object[wgal].irac1 gt 0 then begin
         x(nfluxes) = 35500.
         obs(nfluxes) = object[wgal].irac1
         noise(nfluxes) = 0.05*object[wgal].irac1
         nfluxes = nfluxes+1
      endif
     if object[wgal].irac2 ne 99 and object[wgal].irac2 gt 0 then begin
         x(nfluxes) = 44930.
         obs(nfluxes) = object[wgal].irac2
         noise(nfluxes) = 0.05*object[wgal].irac2
         nfluxes = nfluxes+1
      endif
     if object[wgal].irac3 ne 99 and object[wgal].irac3 gt 0 then begin
         x(nfluxes) = 57310.
         obs(nfluxes) = object[wgal].irac3
         noise(nfluxes) = 0.05*object[wgal].irac3
         nfluxes = nfluxes+1
      endif
     if object[wgal].irac4 ne 99 and object[wgal].irac4 gt 0 then begin
         x(nfluxes) = 78720.
         obs(nfluxes) = object[wgal].irac4
         noise(nfluxes) = 0.05*object[wgal].irac4
         nfluxes = nfluxes+1
      endif

     x = x[0:nfluxes - 1]
     obs = obs[0:nfluxes-1]
     noise = noise[0:nfluxes-1]
  
;     obs = obs*1E-29
;     obs = obs *3E18
;     obs = obs * area
;     obs = obs / 3.827E33
  
;     noise = noise*1E-29
;     noise = noise *3E18
;     noise = noise * area
;     noise = noise / 3.827E33

     ;find the redshifted wavelengths in the model stellar population
     point = fltarr(n_elements(x))
     for i = 0, n_elements(x)-1 do begin
        if i lt 4 then pm = 9. else pm = 90.
        for l = 0, n_elements(lambda) - 1 do begin
           if lambda(l) lt x(i)/(1+z)+ pm and lambda(l) gt x(i)/(1+z) - pm then begin
              point(i) = l
           endif
        endfor
     endfor
     
     ;factor for converting bc03 flux into microjy (as a function of lambda)
     factor = fltarr(n_elements(x)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     factor = factor * (x/(1+z))^2      ;now flux is in microjanskies

     bigfactor = fltarr(n_elements(lambda)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     bigfactor = bigfactor  * (lambda)^2     ;lambda is already shifted, don't need another factor of 1+z
;----------------------------------------------------------------------------
    ;what are the predicted values?

      finalchi = 1E6
      maxmass = 0.1
      minmass = 1E16

      chiarr = fltarr(n_elements(filename)*9)
      massarr = fltarr(n_elements(filename)*9)
      probarr= fltarr(n_elements(filename)*9)
      marr= fltarr(n_elements(filename)*9)
      agearr= fltarr(n_elements(filename)*9)
      count = 0

   ;for each model
      for m = 0, n_elements(filename) - 1 do begin

      ;expected values as a function of age, at the correct wavelengths
      ;for comparison with the observed values
         ex1 = model[0, point ,m]     ;exflux1(point)
         ex3 = model[1,point, m]
         ex5 =  model[2, point ,m]     ;exflux5(point)
         ex7 = model[3,point,m]
         ex10 = model[4, point ,m]     ;exflux10(point)
         ex30 = model[5, point ,m]     ;exflux30(point)
         ex50 = model[6, point ,m]     ;exflux50(point)
         ex80 = model[7, point ,m]     ;exflux80(point)
         ex120 = model[8, point ,m]     ;exflux120(point)

         
         ;convert from solar luminosities per angstrom(bc03 output) to microJy(sensible)          
         ex1 = ex1*factor
         ex3 = ex3*factor
         ex5 = ex5*factor
         ex7 = ex7*factor
         ex10 = ex10*factor
         ex30 = ex30*factor
         ex50 = ex50*factor
         ex80 = ex80*factor
         ex120 = ex120*factor


            ;scale the flux to make the chi-sqrd better?
         start = [5E10]

         result1= MPFITFUN('linear2',ex1,obs, noise, start,/quiet) 
         result3= MPFITFUN('linear2',ex3,obs, noise, start,/quiet) 
         result5= MPFITFUN('linear2',ex5,obs, noise, start,/quiet) 
         result7= MPFITFUN('linear2',ex7,obs, noise, start,/quiet) 
         result10= MPFITFUN('linear2',ex10,obs, noise, start,/quiet) 
         result30= MPFITFUN('linear2',ex30,obs, noise, start,/quiet) 
         result50= MPFITFUN('linear2',ex50,obs, noise, start,/quiet) 
         result80= MPFITFUN('linear2',ex80,obs, noise, start,/quiet) 
         result120= MPFITFUN('linear2',ex120,obs, noise, start,/quiet) 
      ;plot, ob s, result(0)*obs + result(1), thick = 3,psym = 2;, color = colors.green

;         result = [result1, result3,result5,result7,result10,result30,result50,result80,result120]

         


      ;scale the flux by the fitted value
      ;will need this value to scale the stellar mass
         ex1 = result1(0)*ex1
         ex3 = result3(0)*ex3
         ex5 = result5(0)*ex5
         ex7 = result7(0)*ex7
         ex10 = result10(0)*ex10
         ex30 = result30(0)*ex30
         ex50 = result50(0)*ex50
         ex80 = result80(0)*ex80
         ex120 = result120(0)*ex120

 ;        ex = [ex1,ex3,ex5,ex7,ex10,ex30,ex50,ex80,ex120]

      ;determine the chi squared of the fit between two arrays
         df = n_elements(ex1) - 1 ;degrees of freedom

;------
         if maxage ge 1E8 then begin
            residual = (obs - ex1)
            chi = total(residual^(2.0) / noise^(2.0))
            
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result1(0)
            agearr(count) = 0
            marr(count) = m
            count = count + 1
         endif

;------
         if maxage ge 3E8 then begin
            residual = (obs - ex3)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result3(0)
            agearr(count) = 1
            marr(count) = m
            count = count + 1
         endif

;------
         if maxage ge 5E8 then begin
            residual = (obs - ex5)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result5(0)
            agearr(count) = 2
            marr(count) = m
            count = count + 1
        endif
;------
         if maxage ge 7E8 then begin
            residual = (obs - ex7)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result7(0)
            agearr(count) = 3
            marr(count) = m
            count = count + 1
 
         endif

;------
         if maxage ge 1E9 then begin
            residual = (obs - ex10)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result10(0)
            agearr(count) = 4
            marr(count) = m
            count = count + 1
          endif
         
;------
         if maxage ge 3E9 then begin
            residual = (obs - ex30)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result30(0)
            agearr(count) = 5
            marr(count) = m
            count = count + 1
        endif


;------
         if maxage ge 5E9 then begin
            residual = (obs - ex50)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result50(0)
            agearr(count) = 6
            marr(count) = m
            count = count + 1
         endif

;------
         if maxage ge 8E9 then begin
            residual = (obs - ex80)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result80(0)
            agearr(count) = 7
            marr(count) = m
            count = count + 1
        endif

;------
         if maxage ge 12E9 then begin
            residual = (obs - ex120)
            chi = total(residual^(2.0) / noise^(2.0))
            chiarr(count) = chi
            probarr(count) = 1 - chisqr_pdf(chi, df)
            massarr(count) = result120(0)
            agearr(count) = 8
            marr(count) = m
            count = count + 1
          endif




      endfor

      chiarr = chiarr[0:count-1]
      probarr =probarr[0:count-1]
      massarr =massarr[0:count-1]
      agearr = agearr[0:count-1]
      marr = marr[0:count-1]
      
      index = sort(chiarr)
      sortedchi = chiarr[index]
      sortedprob = probarr[index]
      sortedmass = massarr[index]
      sortedage = agearr[index]
      sortedm = marr[index]

      ;determine range of masses based on pm 2% in chisquared

      for mr = 1, 20 do begin
         minlevel = 0.5
         level = 0.02*sortedchi(0)
         if level lt minlevel then level =minlevel
         if sortedchi(mr) le sortedchi(0) +level then begin
            if sortedmass(mr) le sortedmass(0) and sortedmass(mr) lt minmass then minmass = sortedmass(mr)
            if sortedmass(mr) ge sortedmass(0) and sortedmass(mr) gt maxmass then maxmass = sortedmass(mr)
         endif
       endfor
       ;else need to put something in maxmass and minmass
         
      IF maxmass lt 1 and minmass lt 1E16 then maxmass = sortedmass(0) + (sortedmass(0) - minmass)
      if minmass eq 1E16 and maxmass gt 1 then minmass = sortedmass(0) - (maxmass - sortedmass(0))
      if maxmass eq 0.1 and minmass eq 1E16 then begin
         maxmass = sortedmass(0)*2
         minmass = sortedmass(0)*2
      endif

      print, "smallest chi", sortedchi(0), sortedprob(0), sortedmass(0), age(sortedage(0))
      print, "min,max", minmass, maxmass

      plot, x, obs, psym = 2, thick = 3,title = wgal, color = colors.red,/xlog,/ylog
      oplot, lambda*(1+z), model[sortedage(0),*,sortedm(0)]*bigfactor*sortedmass(0) , color = colors.black
      errplot, x,obs -noise, obs +noise, thick = 3


      object[wgal].mass = sortedmass(0)
      object[wgal].masschi = sortedchi(0)
      object[wgal].massprob = sortedprob(0)
      object[wgal].massage = age(sortedage(0))
      object[wgal].model = filename(sortedm(0))
      object[wgal].plusmasserr = maxmass - sortedmass(0)
      object[wgal].minusmasserr = sortedmass(0) - minmass

      endelse
   endif
endfor


print,"Finished at "+systime()

ps_close, /noprint, /noid

save, object, filename='/Users/jkrick/idlbin/object.sav'
end


