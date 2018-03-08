pro stellarmass_plume

  ;device, true=24
  ;device, decomposed=0
  !p.multi = [0, 0, 1]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

ps_start, filename = "/Users/jkrick/Virgo/irac/stellarmass_sed_A.ps"
  
;restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'        ;object
  restore,'/Users/jkrick/bc03/mymodels/bc03_2.sav' ;model
  restore,'/Users/jkrick/bc03/mymodels/lambda_2.sav' ;lambda
  
  restore, '/Users/jkrick/maraston_models/model.2.sav' ;mmodel
  restore, '/Users/jkrick/maraston_models/lambda.sav'  ;mlambda
  restore, '/Users/jkrick/maraston_models/age.sav'     ;age
   
  
  print,"started at "+systime()
  
  readcol,'/Users/jkrick/bc03/mymodels/filelist_2',filename,format="A",/silent
  filename = strcompress(filename + ",0.3,0.5,0.7,1,3,5,8,12")
  
  nlambda = 1220.D
  
;age = fltarr(7)
  bc03age = [0.1,0.3,0.5,0.7,1.,3.,5.,8.,12.]
  
  
;for this plume

  z = .0036
  maxage = galage(z,1000,/silent)

  lumindist = 16.6              ;in Mpc
  lumindist = lumindist * 3.08D24 ;in cm
  area = 0.D
  area = 4*!PI*(lumindist)^2    ;in cm^2
  print, 'lumindist, area', lumindist, area

  x = fltarr(4)                 ;4 wavelengths
  x = [ 4424,5504,35500,44930] ;angstroms

  obs = fltarr(4)
  noise = fltarr(4)
 ; obs = fltarr(2)
 ; noise = fltarr(2)

  nfluxes = 4

;plume A
  a_length = 668.               ;arcsec
  ch1signal = .00067            ;Mjy/sr
  ch1noise = .00036             ;Mjy/sr
  
  ch2signal = .00065            ;Mjy/sr
  ch2noise =.00052              ;Mjy/sr
  
  b_sb = 29.5                   ; mag/sq. arcsec
  b_noisesb = 29.0              ; mag/sq. arcsec
  
  v_sb = 28.6                   ; mag/sq. arcsec
  v_noisesb = 28.5              ; mag/sq. arcsec
  

;plume B
;  a_length = 410.               ;arcsec
;  ch1signal = .00045            ;Mjy/sr
;  ch1noise = .00045             ;Mjy/sr

;  ch2signal = .00055            ;Mjy/sr
;  ch2noise =.00055              ;Mjy/sr

;  b_sb = 30.2                   ; mag/sq. arcsec
;  b_noisesb = 29.0              ; mag/sq. arcsec

;  v_sb = 29.2                  ; mag/sq. arcsec
;  v_noisesb = 28.5              ; mag/sq. arcsec

 ;put everything into microjanskies                              

;convert irac SB to vega mags over the area in question
;ch1-----------------------------------
ch1signal = ch1signal / 3282.8 ; Mjy./sq. degree
ch1signal = ch1signal / ((3.6D3)^2) ; Mjy / sq. arcsec
ch1signal = ch1signal * (a_length^2) ; Mjy
ch1signal = ch1signal * 1E12  ; microjy

ch1noise = ( ch1noise) / 3282.8 ; Mjy./sq. degree
ch1noise = ch1noise / ((3.6D3)^2) ; Mjy / sq. arcsec
ch1noise = ch1noise * (a_length^2) ; Mjy
ch1noise = ch1noise * 1E12  ; microjy

;ch2-----------------------------------
ch2signal = ch2signal / 3282.8 ; Mjy./sq. degree
ch2signal = ch2signal / ((3.6D3)^2) ; Mjy / sq. arcsec
ch2signal = ch2signal * (a_length^2) ; Mjy
ch2signal = ch2signal * 1E12  ; microjy

ch2noise = ( ch2noise) / 3282.8 ; Mjy./sq. degree
ch2noise = ch2noise / ((3.6D3)^2) ; Mjy / sq. arcsec
ch2noise = ch2noise * (a_length^2) ; Mjy
ch2noise = ch2noise * 1E12 ; microjy

;B---------------------------------------
b_vegamag = b_sb -2.5*alog10(a_length^2)
b_abmag = b_vegamag - 0.163
bsignal = magab_to_flux(b_abmag)  ; in erg/s/cm2/Hz
bsignal = (bsignal / 1E-23 ) * 1E6; in microjy

;b_noise_vega = b_noisesb - 2.5*alog10(a_length^2)
;b_noise_ab = b_noise_vega - 0.163
;bnoise = magab_to_flux(b_noise_ab) ; in erg/s/cm2/hz
;bnoise = (bnoise / 1E-23) * 1E6 ; in microjy

;assume 0.5 mag error both pos. and neg.
berr = 0.5
bnoise_pos = (bsignal * 10^(berr/2.5)) - bsignal 
;bnoise_neg = bsignal - (bsignal * 10^(berr/-2.5))

;damn, can't have uneven error bars
;assume the larger one for now

;V---------------------------------------
v_vegamag = v_sb -2.5*alog10(a_length^2)
v_abmag = v_vegamag - 0.044
vsignal = magab_to_flux(v_abmag)  ; in erg/s/cm2/Hz
vsignal = (vsignal / 1E-23 ) * 1E6; in microjy

v_noise_vega = v_noisesb - 2.5*alog10(a_length^2)
v_noise_ab = v_noise_vega - 0.044
vnoise = magab_to_flux(v_noise_ab) ; in erg/s/cm2/hz
vnoise = (vnoise / 1E-23) * 1E6 ; in microjy

;assume 0.5 mag error both pos. and neg.
verr = 0.5
vnoise_pos = (vsignal * 10^(verr/2.5)) - vsignal 
;vnoise_neg = vsignal - (vsignal * 10^(verr/-2.5))

;damn, can't have uneven error bars
;assume the larger one for now

obs = [bsignal, vsignal,ch1signal,ch2signal]
noise = [bnoise_pos, vnoise_pos, ch1noise, ch2noise]
;obs = [bsignal, vsignal]
;noise = [bnoise, vnoise]

print, 'obs', obs
print, 'noise', noise
;----------------------------------------------------------------------------
;bc03                                ;find the redshifted wavelengths in the model stellar population
     point = fltarr(n_elements(x))
     for i = 0, n_elements(x)-1 do begin
        if i lt 2 then pm = 9. else pm = 90.
        for l = 0, n_elements(lambda) - 1 do begin
           ;print, 'bc03', lambda(l),  x(i)/(1+z)+ pm, x(i)/(1+z) - pm
           if lambda(l) lt x(i)/(1+z)+ pm and lambda(l) gt x(i)/(1+z) - pm then begin
              point(i) = l
           endif
        endfor
     endfor
     print, 'point', point
;now the point array holds the index to the bc03 lambda array at the
;right wavelengths for this set of photometry I have given it.


;maraston
     mpoint = fltarr(n_elements(x))
     for i = 0, n_elements(x)-1 do begin
        if x(i) lt 8000 then pm = 10
        if x(i) gt 8000 and x(i) lt 13000 then pm = 25 
        if x(i) gt 13000 and x(i) lt 58000 then pm = 100
        if x(i) gt 58000 then pm = 200
;        print, 'x(i), pm', x(i), pm
        for l = 0, nlambda - 1 do begin
;           print, 'l', l
;           print,'mlambda(l)', mlambda(l), x(i)/(1+z)+ pm,x(i)/(1+z) - pm
           if mlambda(l) lt x(i)/(1+z)+ pm and mlambda(l) gt x(i)/(1+z) - pm then begin
              mpoint(i) = l
           endif
        endfor
     endfor
     print, 'mpoint', mpoint
;----------------------------------------------------------------------------     
;put the models into the same units as the observations
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
      for m = 0, 6 - 1 do begin; n_elements(filename) - 1 do begin
         print, 'working on ', filename(m)
      ;expected values as a function of age, at the correct wavelengths
      ;for comparison with the observed values
         
         for bcnj =0,9-1 do begin  ; for each age
            bcnarr[*,bcnj] = model[bcnj, point ,m] *factor
         
            ;scale the flux to make the chi-sqrd better?
            result(bcnj) = MPFITFUN('linear2',bcnarr[*,bcnj],obs, noise, start,/quiet) ;
 
                                ;scale the flux by the fitted value
                                ;will need this value to scale the stellar mass
 
            bcnarr[*,bcnj] = result(bcnj)*bcnarr[*,bcnj]
            
            print, 'this is the model flux values at our lambda ',  bcnarr[*,bcnj]
            print, 'this is the mass for that model', result(bcnj)

     ;determine the chi squared of the fit between two arrays
            df = n_elements(point) - 1 ;degrees of freedom


            if bcnarr[0,bcnj] ne 0 and maxage ge bc03age(bcnj) then begin

               residual = (obs - bcnarr[*,bcnj])
               chi = total(residual^(2.0) / noise^(2.0))
               print, 'chi', chi
               print, 'prob', 1 - chisqr_pdf(chi, df)
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

      ;print, 'massarr', massarr

      index = sort(chiarr)
      sortedchi = chiarr[index]
      sortedprob = probarr[index]
      sortedmass = massarr[index]
      sortedage = agearr[index]
      sortedm = marr[index]

      ;determine range of masses based on pm 2% in chisquared

      for mr = 1, 6 do begin; 150 do begin
         minlevel = 0.5
         level = 0.2*sortedchi(0)
       ;  if level lt minlevel then level =minlevel
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

      massbc03= sortedmass(0)

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

      ;print, 'msortedchi', msortedchi

      ;determine range of masses based on pm 2% in chisquared

      for mr = 1, 150 do begin
         minlevel = 0.5
         level = 0.2*msortedchi(0)
 ;        if level lt minlevel then level =minlevel
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



      print, "maraston smallest chi", msortedchi(0), msortedprob(0), msortedmass(0);, age(sortedage(0))
      print, "maraston min,max", minmass, maxmass


      massmaraston = msortedmass(0)

;----------------------------------------------------------------------------


      plot, x, obs, psym = 2, thick = 3,/xlog,/ylog,yrange = [100,100000],$
            xtitle = "Angstroms", ytitle = "Flux(microjansky)", charthick=3,xthick=3,ythick=3
      ;errplot, x,obs -noise, obs +noise, thick = 3
      oploterror, x[0:1],obs[0:1],noise[0:1], errthick = 3, psym = 2
      arrow, x[2], obs[2], x[2], obs[2] - noise[2], /data, thick = 3, hthick = 3
      arrow, x[3], obs[3], x[3], obs[3] - noise[3], /data, thick = 3, hthick = 3
      oploterror, x[2:3],obs[2:3],noise[2:3], errthick = 3, psym = 2,/hibar

;bc03
      oplot, lambda*(1+z), model[sortedage(0),*,sortedm(0)]*bigfactor*sortedmass(0), thick = 3

;maraston
      oplot, mlambda*(1+z), mmodel[msortedage(0),*,msortedm(0)]/3.827E33*mbigfactor*msortedmass(0) , linestyle =1, thick = 3

      ;xyouts, 1.1E3, 70, strmid(filename(sortedm(0)), 53,23), charthick = 2
      print, strmid(filename(sortedm(0)), 65,29)


      print,"Finished at "+systime()

      ps_end

end


