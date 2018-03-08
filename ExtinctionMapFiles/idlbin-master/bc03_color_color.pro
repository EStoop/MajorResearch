pro bc03_color_color
ps_start, filename= '/Users/jkrick/Virgo/color_color.ps'

redcolor = FSC_COLOR("Red", !D.Table_Size-2);
orangecolor = FSC_COLOR("Orange", !D.Table_Size-3);
greencolor = FSC_COLOR("Green", !D.Table_Size-4);
bluecolor = FSC_COLOR("Blue", !D.Table_Size-5)
purplecolor = FSC_COLOR("Purple", !D.Table_Size-6)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-7)

filename = '/Users/jkrick/bc03/models/Padova1994/chabrier/bc2003_hr_m52_chab_ssp.color_F124_F129'
readcol, filename,z    ,  LTT  ,   tz   ,   dm     ,    C_rf  ,   C_ne  ,   C_ev  ,    C_AB_rf , C_AB_ne_1_2  ,C_AB_ev_1_2 ,     e_kcor    ,k_cor,format="A",/silent
;plot, tz, C_AB_ne, charthick = 1;, /ylog,/xlog, thick = 1, xrange = [3E4, 6E4], xstyle = 1, yrange = [.1,10.]


filename = '/Users/jkrick/bc03/models/Padova1994/chabrier/bc2003_hr_m52_chab_ssp.color_F130_F131'
readcol, filename,z    ,  LTT  ,   tz   ,   dm     ,    C_rf  ,   C_ne  ,   C_ev  ,    C_AB_rf , C_AB_ne_B_V  ,C_AB_ev_B_V ,     e_kcor    ,k_cor,format="A",/silent


plot, C_AB_ne_B_V, C_AB_ne_1_2, charthick = 1, xtitle = 'B-V (AB)', ytitle = '[3.6-4.5] (AB)', yrange = [-1.0,2.0], title = 'BC03 models, no dust'
xyouts,  C_AB_ne_B_V[130], C_AB_ne_1_2[130], string(tz[130])+ ' Gyr'
xyouts,  C_AB_ne_B_V[100], C_AB_ne_1_2[100], string(tz[100]) + ' Gyr'
xyouts,  C_AB_ne_B_V[50], C_AB_ne_1_2[50], string(tz[50]) + ' Gyr'
xyouts,  C_AB_ne_B_V[1], C_AB_ne_1_2[1], string(tz[1]) + ' Gyr'


;inner, outer
m87_B_V = [1.0, 0.825] - 0.12
m87_1_2 = [-0.56, -0.16]
err_B_V = [ 0.03, 0.075]
err_1_2 = [0.1 , 0.56]

;oplot, m87_B_V, m87_1_2, linestyle = 2
oploterror, m87_B_V, m87_1_2, err_B_V, err_1_2, linestyle = 2
xyouts, m87_B_V[0], m87_1_2[0], 'inner'
xyouts, m87_B_V[1], m87_1_2[1], 'outer'


;----------------------------------------------------------------------------
;why does this look so funny?  try with the bc03 models I have saved and were used in the mass cmd

  restore,'/Users/jkrick/bc03/mymodels/bc03_2.sav' ;model
  restore,'/Users/jkrick/bc03/mymodels/lambda_2.sav' ;lambda
  
  restore, '/Users/jkrick/maraston_models/model.2.sav' ;mmodel
  restore, '/Users/jkrick/maraston_models/lambda.sav'  ;mlambda
  restore, '/Users/jkrick/maraston_models/age.sav'     ;age
   
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

  x = fltarr(4)                 ;3 wavelengths
  x = [ 4424,5504,35500, 45000] ;angstroms

  obs = fltarr(4)
  noise = fltarr(4)
  nfluxes = 4;4


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



;----------------------------------------------------------------------------     
;put the models into the same units as the observations
     ;factor for converting bc03 flux into microjy (as a function of lambda)
     factor = fltarr(n_elements(x)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   

     factor = factor * (x/(1+z))^2      ;now flux is in microjanskies

;bc03
     bigfactor = fltarr(n_elements(lambda)) + 3.827E33/ area / 3E18 / 1E-23 * 1E6   
     bigfactor = bigfactor  * (lambda)^2     ;lambda is already shifted, don't need another factor of 1+z

;----------------------------------------------------------------------------
;bc03
    ;what are the predicted values?

    bcnarr = fltarr(n_elements(point), 9)
    count = 0
    massarr = fltarr(100*9)
    barr = fltarr(100*9)
    varr = fltarr(100*9)
    ch1arr = fltarr(100*9)
    ch2arr = fltarr(100*9)

    m =4 ;m52
    print, 'working on ', filename(m)
    
;    plot, findgen(10), findgen(10), xrange = [13, 17], yrange = [0.5, 1.1], ystyle = 1,/nodata, xtitle = '3.6 (AB)', ytitle = 'B - V (AB)', thick = 3, charthick = 3,xthick = 3, ythick = 3
    
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
          ch2arr(count) =  bcnarr[3,bcnj]*mass
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
       ch2arr =ch2arr[0:count-1]

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
 
       ch2arr = (ch2arr /1E6) *1E-23    ;in erg/s/cm2/hz
       ch2arr = flux_to_magab(ch2arr)   ;in  magab
      
 ;      print, 'massarr', massarr
;       print, 'barr -varr', barr -varr
;       print, 'ch1arr', ch1arr
       
       for a = start, 4E8, size do begin
          b = where(massarr eq a)
 ;         print, 'mass', massarr(b)
 ;         print, 'ch1arr(b)', ch1arr(b)
          oplot, barr(b)-varr(b), ch1arr(b) - ch2arr(b), thick = 3, linestyle = 2
          ;if a gt 2E8 then xyouts, ch1arr(b(8)), barr(b(8)) - varr(b(8)),   string(a, format = '(E10.1)'), orientation = 45, charthick = 3

       endfor
       
     
;       xyouts, 14.4, 0.65, '3Gyr', charthick =3
;       xyouts, 14.9, 0.68, '5Gyr', charthick =3
;       xyouts, 15.15, 0.75, '8Gyr', charthick =3
 ;      xyouts, 15.3, 0.79, '12Gyr', charthick =3
       
       ps_end
    end
