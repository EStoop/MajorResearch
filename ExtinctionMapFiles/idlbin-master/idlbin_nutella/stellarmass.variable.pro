pro stellarmass

device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

!p.multi = [0, 0, 1]
;ps_open, file = "/Users/jkrick/ZPHOT/stellarmass.compare.swire.ps", /portrait, xsize = 4, ysize = 4,/color
ps_open, file = "/Users/jkrick/ZPHOT/variable.ps", /portrait, xsize = 4, ysize = 4,/color

;restore, '/Users/jkrick/idlbin/object.sav'        ;object
restore, '/Users/jkrick/idlbin/variable.sav'        ;object

restore,'/Users/jkrick/bin/bc03/bc03.sav'      ;model
restore,'/Users/jkrick/bin/bc03/lambda.sav'  ;lambda

  restore, '/Users/jkrick/maraston/model.2.sav'    ;mmodel
  restore, '/Users/jkrick/maraston/lambda.sav'   ;mlambda
  restore, '/Users/jkrick/maraston/age.sav'   ;age


print,"started at "+systime()

readcol,'/Users/jkrick/bin/bc03/filelist',filename,format="A",/silent
filename = strcompress(filename + ",0.3,0.5,0.7,1,3,5,8,12")

nlambda = 1220.D

;age = fltarr(7)
bc03age = [0.1,0.3,0.5,0.7,1.,3.,5.,8.,12.]

massmaraston = fltarr(n_elements(variable.rmaga))
massbc03 = fltarr(n_elements(variable.rmaga))
ngals = 0
;for this galaxy
for wgal = 0, n_elements(variable.rmaga) - 1 do begin
   if variable[wgal].prob ge 0. then begin
      ;don't bother calculating stellar mass without a good z calculation
   
      if variable[wgal].gfwhm lt 0.01 and variable[wgal].gfwhm gt 0 then begin
         ;don't calculate for stars
         print, "rejecting a star"
      endif else begin
      print, "galaxy ", wgal, variable[wgal].zphot

      z = variable[wgal].zphot
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
      if variable[wgal].umaga lt 30 and variable[wgal].umaga gt 10 then begin
         x(nfluxes) = 3540.
         obs(nfluxes) = 1E6*(10^((variable[wgal].umaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((variable[wgal].umaga - 8.926)/(-2.5))) - (10^((variable[wgal].umaga +variable[wgal].umagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes + 1
      endif
      if variable[wgal].gmaga lt 30 and variable[wgal].gmaga gt 10 then begin
         x(nfluxes) = 4660.
         obs(nfluxes) = 1E6*(10^((variable[wgal].gmaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((variable[wgal].gmaga - 8.926)/(-2.5))) - (10^((variable[wgal].gmaga +variable[wgal].gmagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes + 1
      endif
      if variable[wgal].rmaga lt 30 and variable[wgal].rmaga gt 10 then begin
         x(nfluxes) = 6255.
         obs(nfluxes) =1E6*(10^((variable[wgal].rmaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((variable[wgal].rmaga - 8.926)/(-2.5))) - (10^((variable[wgal].rmaga +variable[wgal].rmagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes + 1
      endif
      if variable[wgal].imaga lt 30 and variable[wgal].imaga gt 10 then begin
         x(nfluxes) = 7680.
         obs(nfluxes) =1E6*(10^((variable[wgal].imaga - 8.926)/(-2.5)))
         noise(nfluxes) = 1E6*((10^((variable[wgal].imaga - 8.926)/(-2.5))) - (10^((variable[wgal].imaga +variable[wgal].imagerra- 8.926)/(-2.5))))
         nfluxes = nfluxes+1
      endif
      if variable[wgal].flamjmag lt 30 and variable[wgal].flamjmag gt 10 then begin
         x(nfluxes) = 12490.
         obs(nfluxes) =1E6*1594.*10^(variable[wgal].flamjmag/(-2.5))
         noise(nfluxes) = 1E6*1594.*(10^(variable[wgal].flamjmag/(-2.5)) - (10^((variable[wgal].flamjmag+ variable[wgal].flamjmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif
      if variable[wgal].wircjmag lt 30 and variable[wgal].wircjmag gt 10 then begin
         x(nfluxes) = 12500
         obs(nfluxes) =1E6*1594.*10^(variable[wgal].wircjmag/(-2.5))
         noise(nfluxes) = 1E6*1594.*(10^(variable[wgal].wircjmag/(-2.5)) - (10^((variable[wgal].wircjmag+ variable[wgal].wircjmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif

     if variable[wgal].tmassj ne 99 and variable[wgal].tmassj gt 0 then begin
         x(nfluxes) = 12400.
         obs(nfluxes) = 1E6*variable[wgal].tmassj
         noise(nfluxes) = 0.05*1E6*variable[wgal].tmassj
         nfluxes = nfluxes+1
      endif
      if variable[wgal].wirchmag lt 30 and variable[wgal].wirchmag gt 10 then begin
         x(nfluxes) = 16350
         obs(nfluxes) =1E6*1024.*10^(variable[wgal].wirchmag/(-2.5))
         noise(nfluxes) = 1E6*1024.*(10^(variable[wgal].wirchmag/(-2.5)) - (10^((variable[wgal].wirchmag+ variable[wgal].wirchmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif
     if variable[wgal].tmassh ne 99 and variable[wgal].tmassh gt 0 then begin
         x(nfluxes) = 16500.
         obs(nfluxes) = 1E6*variable[wgal].tmassh
         noise(nfluxes) = 0.05*1E6*variable[wgal].tmassh
         nfluxes = nfluxes+1
      endif
      if variable[wgal].wirckmag lt 30 and variable[wgal].wirckmag gt 10 then begin
         x(nfluxes) = 21500
         obs(nfluxes) =1E6*666.8*10^(variable[wgal].wirckmag/(-2.5))
         noise(nfluxes) = 1E6*666.8*(10^(variable[wgal].wirckmag/(-2.5)) - (10^((variable[wgal].wirckmag+ variable[wgal].wirckmagerr)/(-2.5))))
         nfluxes = nfluxes+1
      endif
     if variable[wgal].tmassk ne 99 and variable[wgal].tmassk gt 0 then begin
         x(nfluxes) = 21600.
         obs(nfluxes) = 1E6*variable[wgal].tmassk
         noise(nfluxes) = 0.05*1E6*variable[wgal].tmassk
         nfluxes = nfluxes+1
      endif
      if variable[wgal].irac1 ne 99 and variable[wgal].irac1 gt 0 then begin
         x(nfluxes) = 35500.
         obs(nfluxes) = variable[wgal].irac1
         noise(nfluxes) = 0.05*variable[wgal].irac1
         nfluxes = nfluxes+1
      endif
     if variable[wgal].irac2 ne 99 and variable[wgal].irac2 gt 0 then begin
         x(nfluxes) = 44930.
         obs(nfluxes) = variable[wgal].irac2
         noise(nfluxes) = 0.05*variable[wgal].irac2
         nfluxes = nfluxes+1
      endif
     if variable[wgal].irac3 ne 99 and variable[wgal].irac3 gt 0 then begin
         x(nfluxes) = 57310.
         obs(nfluxes) = variable[wgal].irac3
         noise(nfluxes) = 0.05*variable[wgal].irac3
         nfluxes = nfluxes+1
      endif
     if variable[wgal].irac4 ne 99 and variable[wgal].irac4 gt 0 then begin
         x(nfluxes) = 78720.
         obs(nfluxes) = variable[wgal].irac4
         noise(nfluxes) = 0.05*variable[wgal].irac4
         nfluxes = nfluxes+1
      endif
 
   if nfluxes gt 3 then begin   ;don't bother with only 3 points

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

 

;----------------------------------------------------------------------------
;bc03                                ;find the redshifted wavelengths in the model stellar population
     point = fltarr(n_elements(x))
     for i = 0, n_elements(x)-1 do begin
        if i lt 4 then pm = 9. else pm = 90.
        for l = 0, n_elements(lambda) - 1 do begin
           if lambda(l) lt x(i)/(1+z)+ pm and lambda(l) gt x(i)/(1+z) - pm then begin
              point(i) = l
           endif
        endfor
     endfor
;----------------------------------------------------------------------------     
;maraston
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
;----------------------------------------------------------------------------     

     ;factor for converting bc03 flux into microjy (as a function of lambda)
     factor = fltarr(n_elements(x)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     factor = factor * (x/(1+z))^2      ;now flux is in microjanskies

;bc03
     bigfactor = fltarr(n_elements(lambda)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     bigfactor = bigfactor  * (lambda)^2     ;lambda is already shifted, don't need another factor of 1+z
;maraston
    mbigfactor = fltarr(n_elements(mlambda)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
    mbigfactor = mbigfactor  * (mlambda)^2 ;lambda is already shifted, don't need another factor of 1+z

;----------------------------------------------------------------------------
;bc03
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
      bcnarr = fltarr(n_elements(point), 9)
      start = [6E10]
      result = fltarr(9)

   ;for each model
      for m = 0, n_elements(filename) - 1 do begin

      ;expected values as a function of age, at the correct wavelengths
      ;for comparison with the observed values
         
         for bcnj =0,9-1 do begin
            bcnarr[*,bcnj] = model[bcnj, point ,m] *factor
 
            ;scale the flux to make the chi-sqrd better?
            result(bcnj) = MPFITFUN('linear2',bcnarr[*,bcnj],obs, noise, start,/quiet) 
 
                                ;scale the flux by the fitted value
                                ;will need this value to scale the stellar mass
 
            bcnarr[*,bcnj] = result(bcnj)*bcnarr[*,bcnj]


     ;determine the chi squared of the fit between two arrays
         df = n_elements(point) - 1 ;degrees of freedom


            if bcnarr[0,bcnj] ne 0 and maxage ge bc03age(bcnj) then begin

               residual = (obs - bcnarr[*,bcnj])
               chi = total(residual^(2.0) / noise^(2.0))
               
               chiarr(count) = chi
               probarr(count) = 1 - chisqr_pdf(chi, df)
               massarr(count) = result(bcnj)
               agearr(count) = bcnj
               marr(count) = m
               
               count = count + 1
            endif

         endfor

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

      print, "smallest chi", sortedchi(0), sortedprob(0), sortedmass(0), bc03age(sortedage(0))
      print, "min,max", minmass, maxmass

      massbc03(ngals) = sortedmass(0)
;----------------------------------------------------------------------------
;      maraston
    ;what are the predicted values?

      finalchi = 1E6
      maxmass = 0.1
      minmass = 1E16

      mchiarr = fltarr(67*68)
      mmassarr = fltarr(67*68)
      mprobarr= fltarr(67*68)
      mmarr= fltarr(67*68)
      magearr= fltarr(67*68)
      mcount = 0

      narr = fltarr(n_elements(mpoint), 67)
      start = [6E10]
      mresult = fltarr(67)

   ;for each model
      for m = 0, 68 - 1 do begin

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


      massmaraston(ngals) = msortedmass(0)

;----------------------------------------------------------------------------


      plot, x, obs, psym = 2, thick = 3, /xlog,/ylog,$
            title = strcompress(string(variable[wgal].zphot) + string(variable[wgal].prob) + string(minmass) + string(maxmass)), xtitle = "Angstroms", ytitle = "flux(microjansky)", charthick=3,xthick=3,ythick=3
   ;         title=strcompress(string(variable[wgal].ra) +string(variable[wgal].dec) +string(variable[wgal].zphot) + string(variable[wgal].prob)), xtitle = "Angstroms", ytitle = "flux"
    


  errplot, x,obs -noise, obs +noise, thick = 3

;bc03
      oplot, lambda*(1+z), model[sortedage(0),*,sortedm(0)]*bigfactor*sortedmass(0) , color = colors.black

;maraston
      oplot, mlambda*(1+z), mmodel[msortedage(0),*,msortedm(0)]/3.827E33*mbigfactor*msortedmass(0) , color = colors.blue

;      xyouts, 1.1E3, 70, strmid(filename(sortedm(0)), 53,23), charthick = 2

     variable[wgal].mass = sortedmass(0)
     variable[wgal].masschi = sortedchi(0)
     variable[wgal].massprob = sortedprob(0)
     variable[wgal].massage = bc03age(sortedage(0))
     variable[wgal].model = filename(sortedm(0))
     variable[wgal].plusmasserr = maxmass - sortedmass(0)
     variable[wgal].minusmasserr = sortedmass(0) - minmass
      
      ngals = ngals + 1
   endif

      endelse
   endif
endfor


;massbc03 = massbc03(0:ngals - 1)
;massmaraston = massmaraston(0:ngals-1)
;save, massbc03, filename='/Users/jkrick/maraston/massbc03.sav'
;save, massmaraston, filename='/Users/jkrick/maraston/massmaraston.sav'
;
;plot, alog10(massmaraston), alog10(massbc03), xtitle = "maraston mass", $
;      ytitle = "bc03 mass", xrange = [4,15],yrange=[4,15], psym = 2,$
;      xstyle = 1, ystyle = 1
;oplot, findgen(20), findgen(20), thick = 3
print,"Finished at "+systime()

ps_close, /noprint, /noid

save, variable, filename='/Users/jkrick/idlbin/variable.sav'
end


