pro stellarmass_plume

  ;device, true=24
  ;device, decomposed=0
  !p.multi = [0, 0, 1]
  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

ps_start, filename = "/Users/jkrick/Virgo/irac/test_stellarmass_A.ps"
  
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

;  x = fltarr(4)                 ;4 wavelengths
;  x = [ 4424,5504,35500,44930] ;angstroms
  x = fltarr(3)                 ;3 wavelengths
  x = [ 4424,5504,35500] ;angstroms

;  obs = fltarr(4)
;  noise = fltarr(4)
  obs = fltarr(3)
  noise = fltarr(3)

  nfluxes = 3;4

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
b_noise = b_noisesb - 2.5*alog10(a_length^2)
bsignal = magab_to_flux(b_vegamag)  ; know this is AB and not vega, but I just want the delta; in erg/s/cm2/Hz
bsignal = (bsignal / 1E-23 ) * 1E6; in microjy

bnoise = magab_to_flux(b_noise) ; in erg/s/cm2/hz
bnoise = (bnoise / 1E-23) * 1E6 ; in microjy

;V---------------------------------------
v_vegamag = v_sb -2.5*alog10(a_length^2)
v_noise = v_noisesb - 2.5*alog10(a_length^2)
vsignal = magab_to_flux(v_vegamag)  ; know this is AB and not vega, but I just want the delta; in erg/s/cm2/Hz
vsignal = (vsignal / 1E-23 ) * 1E6; in microjy

vnoise = magab_to_flux(v_noise) ; in erg/s/cm2/hz
vnoise = (vnoise / 1E-23) * 1E6 ; in microjy

;obs = [bsignal, vsignal,ch1signal,ch2signal]
;noise = [bnoise, vnoise, ch1noise, ch2noise]

obs = [bsignal, vsignal,ch1signal]
noise = [bnoise, vnoise, ch1noise]

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

       massarr = fltarr(n_elements(filename)*9)
      marr= fltarr(n_elements(filename)*9)
      barr= fltarr(n_elements(filename)*9)
      varr= fltarr(n_elements(filename)*9)
      ch1arr= fltarr(n_elements(filename)*9)
      agearr= fltarr(n_elements(filename)*9)
      count = 0
      bcnarr = fltarr(n_elements(point), 9)
      start = [6E10]
      result = fltarr(9)

      m = 3
      print, 'working on ', filename(m)

          
      bcnj = 6                  ; set to a specific age
      bcnarr[*,bcnj] = model[bcnj, point ,m] *factor
         
      modelarr = bcnarr
      ;change the b flux values and see what that does to the mass
      for bflux = 0, 600, 100 do begin
         obs = [bflux, obs[1], obs[2]]
         print, 'inside obs', obs
         print, 'inside noise', noise
            ;scale the flux to make the chi-sqrd better?
         result = MPFITFUN('linear2',bcnarr[*,bcnj],obs, noise, start) ;
 
                                ;scale the flux by the fitted value
                                ;will need this value to scale the stellar mass
 
         modelarr[*,bcnj] = result*bcnarr[*,bcnj]
            
         print, 'this is the model flux values at our lambda ',  modelarr[*,bcnj] 
         print, 'this is the mass for that model', result

     ;determine the chi squared of the fit between two arrays
     ;       df = n_elements(point) - 1 ;degrees of freedom


         if bcnarr[0,bcnj] ne 0  then begin ;and maxage ge bc03age(bcnj)

            massarr(count) = result
            barr(count) = bcnarr[0,bcnj]
            varr(count) = bcnarr[1,bcnj]
            ch1arr(count) = bcnarr[2,bcnj]
         ;agearr(count) = bcnj
         ;marr(count) = m
            count = count + 1

         endif

      endfor

      massarr =massarr[0:count-1]
      barr =barr[0:count-1]
      varr =varr[0:count-1]
      ch1arr =ch1arr[0:count-1]

 ;     agearr = agearr[0:count-1]
 ;     marr = marr[0:count-1]

      print, 'massarr', massarr

 
;------------

      ps_end

end


