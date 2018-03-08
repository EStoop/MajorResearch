pro maraston_test

close, /all

;device, true=24
;device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/maraston/sed.ps", /portrait, xsize = 4, ysize = 4,/color

;readcol,'/Users/jkrick/maraston/Sed_Mar05_SSP_Salpeter/sed.ssz002.rhb',age,metal,lambda,flux,format="A"
nlambda = 1220.D
;plot, alog10(lambda(0:nlambda)), flux(0:nlambda),/ylog,yrange=[1E23,1E34], ystyle = 1, xrange = [2, 5]

;for i = 1.D, 66, 1 do begin

;   oplot, alog10(lambda((nlambda*i)+(i):nlambda*(i+1) +(i))), flux((nlambda*i)+(i):nlambda*(i+1) +(i));
;endfor




restore, '/Users/jkrick/idlbin/object.sav'
  restore, '/Users/jkrick/maraston/model.sav'
  restore, '/Users/jkrick/maraston/lambda.sav'
  restore, '/Users/jkrick/maraston/age.sav'


print,"started at "+systime()

;for this galaxy
for wgal = 0, 10 do begin;n_elements(object.rmaga) - 1 do begin
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
     mpoint = fltarr(n_elements(x))
     for i = 0, n_elements(x)-1 do begin
        if x(i) lt 8000 then pm = 10
        if x(i) gt 8000 and x(i) lt 13000 then pm = 25 
        if x(i) gt 13000 and x(i) lt 58000 then pm = 50
        if x(i) gt 58000 then pm = 200
        for l = 0, nlambda - 1 do begin
           if mlambda(l) lt x(i)/(1+z)+ pm and mlambda(l) gt x(i)/(1+z) - pm then begin
              mpoint(i) = l
           endif
        endfor
     endfor
     
     ;factor for converting bc03 flux into microjy (as a function of lambda)
     factor = fltarr(n_elements(x)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     factor = factor * (x/(1+z))^2      ;now flux is in microjanskies

     mbigfactor = fltarr(n_elements(mlambda)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     mbigfactor = mbigfactor  * (mlambda)^2     ;lambda is already shifted, don't need another factor of 1+z
;----------------------------------------------------------------------------
    ;what are the predicted values?

      finalchi = 1E6
      maxmass = 0.1
      minmass = 1E16

      mchiarr = fltarr(67*11)
      mmassarr = fltarr(67*11)
      mprobarr= fltarr(67*11)
      mmarr= fltarr(67*11)
      magearr= fltarr(67*11)
      mcount = 0

      narr = fltarr(n_elements(mpoint), 67)
      start = [6E10]
      mresult = fltarr(67)

   ;for each model
      for m = 0, 11 - 1 do begin

      ;expected values as a function of age, at the correct wavelengths
      ;for comparison with the observed values

         for nj = 0,67 - 1 do begin
            narr[*,nj] = mmodel[nj, mpoint ,m] /3.827E33*factor
                                ;convert from erg per s per angstrom(maraston output) to microJy(sensible)  
                                ;keep the same factor from bc03 and convert maraston to bc03 though solar lum.

            ;scale the flux to make the chi-sqrd better?
            mresult(nj) = MPFITFUN('linear2',narr[*,nj],obs, noise, start,/quiet) 
 
                                ;scale the flux by the fitted value
                                ;will need this value to scale the stellar mass


            narr[*,nj] = mresult(nj)*narr[*,nj]


                                ;determine the chi squared of the fit between two arrays
                                ;may want to add in, not to consider the ones older than the universe at that redshift
            df = n_elements(mpoint) - 1 ;degrees of freedom
       
            if narr[0,nj] ne 0 and maxage ge age(nj)*1E9 then begin

               residual = (obs - narr[*,nj])
               chi = total(residual^(2.0) / noise^(2.0))
               
               mchiarr(mcount) = chi
               mprobarr(mcount) = 1 - chisqr_pdf(chi, df)
               mmassarr(mcount) = mresult(nj)
               magearr(mcount) = nj
               mmarr(mcount) = m
               
               mcount = mcount + 1
            endif

         endfor



      endfor

      mchiarr = mchiarr[0:mcount-1]
      mprobarr =mprobarr[0:mcount-1]
      mmassarr =mmassarr[0:mcount-1]
      magearr = magearr[0:mcount-1]
      mmarr = mmarr[0:mcount-1]


      index = sort(mchiarr)
      msortedchi = mchiarr[index]
      msortedprob = mprobarr[index]
      msortedmass = mmassarr[index]
      msortedage = magearr[index]
      msortedm = mmarr[index]



      ;determine range of masses based on pm 2% in chisquared

      for mr = 1, 20 do begin
         minlevel = 0.5
         level = 0.02*msortedchi(0)
         if level lt minlevel then level =minlevel
         if msortedchi(mr) le msortedchi(0) +level then begin
            if msortedmass(mr) le msortedmass(0) and msortedmass(mr) lt minmass then minmass = msortedmass(mr)
            if msortedmass(mr) ge msortedmass(0) and msortedmass(mr) gt maxmass then maxmass = msortedmass(mr)
         endif
      endfor

       ;else need to put something in maxmass and minmass
         
      IF maxmass lt 1 and minmass lt 1E16 then maxmass = msortedmass(0) + (msortedmass(0) - minmass)
      if minmass eq 1E16 and maxmass gt 1 then minmass = msortedmass(0) - (maxmass - msortedmass(0))
      if maxmass eq 0.1 and minmass eq 1E16 then begin
         maxmass = msortedmass(0)*2
         minmass = msortedmass(0)*2
      endif



      print, "smallest chi", msortedchi(0), msortedprob(0), msortedmass(0);, age(sortedage(0))
      print, "min,max", minmass, maxmass

      plot, x, obs, psym = 2, thick = 3,title = wgal, color = colors.red,/xlog,/ylog

      oplot, mlambda*(1+z), mmodel[msortedage(0),*,msortedm(0)]/3.827E33*mbigfactor*msortedmass(0) , color = colors.black
      errplot, x,obs -noise, obs +noise, thick = 3


      object[wgal].mass = msortedmass(0)
      object[wgal].masschi = msortedchi(0)
      object[wgal].massprob = msortedprob(0)
      object[wgal].massage = age(msortedage(0))
;      object[wgal].model = filename(sortedm(0))
      object[wgal].plusmasserr = maxmass - msortedmass(0)
      object[wgal].minusmasserr = msortedmass(0) - minmass

      endelse
   endif
endfor


print,"Finished at "+systime()

ps_close, /noprint, /noid

;save, object, filename='/Users/jkrick/idlbin/object.maraston.sav'
end






