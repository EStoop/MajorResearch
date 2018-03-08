pro print_bc03
;-------------------------------------------------------------------------------------------------------
;This function prints out an input file for python code to run bc03



metalarr =[0.0001,0.0004,0.004,0.008,00.2,0.05]
metalnamearr = ["m22","m32","m42","m52","m62","m72"]
dustarr = ["Y","N"]
dustnamearr = ["d", "nd"]
tauvarr = [0.5,0.6,0.7,0.8,0.9,1.0]
tauvnamearr=[string(05),string(06),string(07),string(08),string(09),string(10)]
sfharr=["ssp","exponential"]
tauarr = [0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1.0,3.0,5.0]
taunamearr=[string(001),string(003),string(005),string(007),string(01),string(03),string(05),string(07),string(10),string(30),string(50)]


for sind = 0 , n_elements(sfharr) -1 do begin
   sfh = sfharr(sind)
   if sind eq 1 then begin
    ;exponential
      for tauind = 0, n_elements(tauarr) - 1 do begin
         tau = tauarr(tauind)
         for mind = 0, n_elements(metalarr) - 1 do begin 
            metal = metalarr(mind)
            for dind = 0, n_elements(dustarr) - 1 do begin
               if dind eq 0 then begin
                  for tauvind = 0, n_elements(tauvarr) - 1 do begin
                     tauv = tauvarr(tauvind)
                     outfile = strcompress("/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_"+$
                                           metalnamearr(mind)+"_salp_"+dustnamearr(dind)+"_"+tauvnamearr(tauvind)+ $
                                           "_"+strmid(sfharr,0,3)+"_"+taunamearr(tauind),/remove_all)
                     print, outfile

                  endfor
               endif else begin
                  tauv = 0
                  tauvind = 0
                  outfile = strcompress("/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_"+$
                                        metalnamearr(mind)+"_salp_"+dustnamearr(dind)+"_"+tauvnamearr(tauvind)+ $
                                        "_"+strmid(sfharr,0,3)+"_"+taunamearr(tauind),/remove_all)
                  print, outfile
               endelse
            endfor
         endfor
      endfor
   endif else begin
      tau = 1.0
      for mind = 0, n_elements(metalarr) - 1 do begin 
         metal = metalarr(mind)
         for dind = 0, n_elements(dustarr) - 1 do begin
            if dind eq 0 then begin
               for tauvind = 0, n_elements(tauvarr) - 1 do begin
                  tauv = tauvarr(tauvind)
                  tauind = 0
                  outfile = strcompress("/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_"+$
                                        metalnamearr(mind)+"_salp_"+dustnamearr(dind)+"_"+tauvnamearr(tauvind)+ $
                                        "_"+strmid(sfharr,0,3)+"_"+taunamearr(tauind),/remove_all)
                  print, outfile

               endfor
            endif else begin
               tauv = 0
               tauvind = 0
               tauind = 0
               outfile = strcompress("/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_"+$
                                     metalnamearr(mind)+"_salp_"+dustnamearr(dind)+"_"+tauvnamearr(tauvind)+ $
                                     "_"+strmid(sfharr,0,3)+"_"+taunamearr(tauind),/remove_all)
               print, outfile
            endelse
         endfor
      endfor
   endelse
endfor


close, /all

END


;                 outfile = strcompress("/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_"+$
;                                        metalnamearr(mind)+"_salp_"+dustnamearr(dind)+"_"+tauvnamearr(tind)+ $
;                                        "_"+strmid(sfharr,0,3)+"_"+taunamearr(tauind),/remove_all)
;                  print, outfile
                 openw, outlun, "/Users/jkrick/bin/bc03/inputfile.auto", /GET_LUN
           
                 printf, outlun, "# Compute stellar evolution of a composite stellar population"
                 printf, outlun, "ROOTDIR		/Users/jkrick/bin/bc03/ " ;#Root BC03 directory
                 printf, outlun, "IMF		salpeter	"                 ;#salpeter/chabrier
                 printf, outlun, "RESOLUTION	lr 	"                         ;	#high = hr, low = lr
                 printf, outlun, "METALLICITY ", metal                           ;	0.05
                 printf, outlun, "TRACKS		Padova1994 "              ;	#Padova1994 / Padova2000
                 printf, outlun, "DUST	", 	 dust                             ;#Y/N. 
                 printf, outlun, "TAU_V 	", tauv                           ; 	#total effective attenuation optical depth"
                 printf, outlun, "MU		0.3	"                         ;#fraction of tau_V arising from the ambient ISM
                 printf, outlun, "SFH	", sfh                                    ;	exponential          # [exponential/SSP/singleburst]
                 printf, outlun, "TAU	", tau                                    ;	0.7	#SFH exponential:tau , SFH singleburst:DeltaT"
                 printf, outlun, "GASRECYCLE	N	"                         ;#Y/N "
                 printf, outlun, "EPSILON		0.001 "                   ;	#recycled fraction of gas (from 0.001 to 1)"
                 printf, outlun, "TCUTSFR		20	"                 ;#Set SFR = 0 at time > TcutSFR (Gyr)"
                 printf, outlun, "CSPOUTFILE", outfile ;	models/Padova1994/salpeter/lr_m72_salp_dust_exp_07	#basename of output files"
                 printf, outlun, "#_____________                      _______________"
                 printf, outlun, "#Extract a spectrum at a given age"
                 printf, outlun, "EXTRACT_AGE	0.1,0.3,0.5,0.7,1,3,5,8,12	#output age(s) (Gyr), max=24 ages; N=do not extract "
                 printf, outlun, "SPECRANGEMIN	500 	#minimum wavelength range (Angstrom)"
                 printf, outlun, "SPECRANGEMAX	81000	#maximum wavelength range (Angstrom)"
                 printf, outlun, "FNU_LAMBDA	lambda	#nu/lambda "
                 printf, outlun, "F0SCALING	N 	#Normalization flux at wavelength W0. N= no scaling"
                 printf, outlun, "W0SCALING	4000	#W0"
                 printf, outlun, "#_____________                      _______________"
                 
                 
                 close, outlun
                 free_lun, outlun
                
                
                 spawn, "python /Users/jkrick/bin/bc03/runBC03reduced.py /Users/jkrick/bin/bc03/inputfile.auto"
                  
