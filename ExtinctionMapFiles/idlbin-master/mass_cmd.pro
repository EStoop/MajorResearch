pro mass_cmd

  ;device, true=24
  ;device, decomposed=0
 ; !p.multi = [0, 0, 1]
  ;redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  ;bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

;ps_start, filename = "/Users/jkrick/Virgo/irac/mass_cmd.ps"
  ps_open, filename='/Users/jkrick/Virgo/irac/mass_cmd.ps',/portrait,/square,/color

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

    bcnarr = fltarr(n_elements(point), 9)
    count = 0
    massarr = fltarr(100*9)
    barr = fltarr(100*9)
    varr = fltarr(100*9)
    ch1arr = fltarr(100*9)
    
    m =3 ;m52
    print, 'working on ', filename(m)
    
    plot, findgen(10), findgen(10), xrange = [13, 17], yrange = [0.5, 1.1], ystyle = 1,/nodata, xtitle = '3.6 (AB)', ytitle = 'B - V (AB)', thick = 3, charthick = 3,xthick = 3, ythick = 3
    
    for bcnj = 0, 8 do begin
       print, 'working on age', bcnj
;      bcnj = 6                  ; set to a specific age
       bcnarr[*,bcnj] = model[bcnj, point ,m] *factor
       modelarr = bcnarr
       
                                ;plot, x, bcnarr, xrange = [1E3, 1E5], yrange = [1E2, 1E5], /xlog, /ylog, /nodata
       start = 1E8
       size = 1E8
       for mass = start, 5E9, size do begin
         ; print, 'working on mass', mass
          modelarr[*,bcnj] = mass*bcnarr[*,bcnj]
          
;            print, 'this is the model flux values at our lambda ',  modelarr[*,bcnj] 
;            print, 'this is the mass for that model', mass
;          print, 'count', count, n_elements(massarr)
          massarr(count) = mass
          barr(count) = bcnarr[0,bcnj]*mass
          varr(count) = bcnarr[1,bcnj]*mass
          ch1arr(count) = bcnarr[2,bcnj]*mass
                                ;agearr(count) = bcnj
                                ;marr(count) = m
          count = count + 1
          
                                ;oplot, x, bcnarr[*,bcnj]*mass
          
       endfor
     endfor
      
       massarr =massarr[0:count-1]
       barr =barr[0:count-1]
       varr =varr[0:count-1]
       ch1arr =ch1arr[0:count-1]
       
                                ;     agearr = agearr[0:count-1]
                                ;     marr = marr[0:count-1]
;------------
       
;need to convert flux back to mags so I can make a color mag. plot.
       barr = (barr /1E6) *1E-23    ;in erg/s/cm2/hz
       barr = flux_to_magab(barr)   ;in  magab
       
       varr = (varr /1E6) *1E-23    ;in erg/s/cm2/hz
       varr = flux_to_magab(varr)   ;in  magab
       
       ch1arr = (ch1arr /1E6) *1E-23    ;in erg/s/cm2/hz
       ch1arr = flux_to_magab(ch1arr)   ;in  magab
       
 ;      print, 'massarr', massarr
;       print, 'barr -varr', barr -varr
;       print, 'ch1arr', ch1arr
       
       for a = start, 4E8, size do begin
          b = where(massarr eq a)
 ;         print, 'mass', massarr(b)
 ;         print, 'ch1arr(b)', ch1arr(b)
          oplot, ch1arr(b), barr(b)-varr(b), thick = 3
          if a gt 2E8 then xyouts, ch1arr(b(8)), barr(b(8)) - varr(b(8)),   string(a, format = '(E10.1)'), orientation = 45, charthick = 3

       endfor
       
           for a = 5E8, 4.5E9,5E8 do begin
              b = where(massarr eq a)
;              print, 'mass', massarr(b)
;              print, 'ch1arr(b)', ch1arr(b)
              oplot, ch1arr(b), barr(b)-varr(b), thick = 3
              
;          modelname =  strmid(filename(m), 65,19) + string(bc03age(bcnj))
              xyouts, ch1arr(b(8)), barr(b(8)) - varr(b(8)),   string(a, format = '(E10.1)'), orientation = 45, charthick = 3
;          if a eq 1.0E9 then begin
;             xyouts, ch1arr(b) - 0.3, barr(b) - varr(b) - 0.01, string(fix(bc03age)) + 'Gyr', charthick =3
;          endif
              
           endfor
           
           xyouts, 14.4, 0.65, '3Gyr', charthick =3
           xyouts, 14.9, 0.68, '5Gyr', charthick =3
           xyouts, 15.15, 0.75, '8Gyr', charthick =3
           xyouts, 15.3, 0.79, '12Gyr', charthick =3
           
;now for plume A itself
microjy = 7027
erg = (microjy /1E6) *1E-23    ;in erg/s/cm2/hz
x_ch1 = flux_to_magab(erg)   ;in  magab
xerr = 0.65
y_bv = 0.9-0.12 ; include conversion from Vega to AB
yerr = 0.15

oploterror, x_ch1,y_bv, xerr, yerr, errthick = 3, psym = 2, thick = 3, /lobar
oploterror, x_ch1,y_bv, yerr, errthick = 3, psym = 2, thick = 3
arrow,x_ch1,   y_bv,x_ch1 + xerr,y_bv, /data, thick = 3, hthick = 3
print, x_ch1, 'x_ch1',  x_ch1 + xerr, y_bv

xyouts, 14.35, 0.92, 'A', charthick = 3

;and plume B

microjy = 1778
erg = (microjy /1E6) *1E-23    ;in erg/s/cm2/hz
x_ch1 = flux_to_magab(erg)   ;in  magab
print, 'x', x_ch1
xerr = 0.8
y_bv = 1.0-0.12 ; include conversion from Vega to AB
yerr = 0.2
oploterror, x_ch1,y_bv, xerr, yerr, errthick = 3, psym = 2, thick = 3,/lobar
oploterror, x_ch1,y_bv, yerr, errthick = 3, psym = 2, thick = 3
arrow,x_ch1,   y_bv,x_ch1 + xerr,y_bv, /data, thick = 3, hthick = 3

xyouts, 15.85, 1.02, 'B', charthick = 3

;draw an arrow that will indicate the direction of a change in metallicity.
arrow,13.66, 0.8, 13.62,0.9, /data, thick = 3, hthick = 3
xyouts, 13.3, 0.9, 'Increase Metallicity', charthick = 3

;ps_end
ps_close, /noprint,/noid

end


