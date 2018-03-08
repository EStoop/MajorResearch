  

function findobject, rawant, decwant
;device, true=24
;device, decomposed=0

close, /all
;rawant = dblarr(1)
;decwant  = rawant
;rawant = [264.96030D,264.86871D]
;decwant=[68.87715D,68.886864D]


;rawant = [265.41046,265.41083,265.42557,265.39478,265.38852,264.58295,264.58136]
     
;decwant=[68.831825,68.831970,68.837524,68.854172,68.859879,68.869736,68.885330]


;matching 70 micron by hand to acs

;rawant = [265.0260925,264.9787598,265.3198242,265.2285767,265.1435852,265.1298218,265.0886536,264.9353638,264.8025818,264.8103638,264.7998962,264.7887878,264.8887329,265.0636597,265.1195984,264.9949951,264.8736267,264.9096985,264.6955566,264.6903381,265.3695984]

;decwant = [69.1017456,69.0819473,69.0449753,69.0356827,69.0588303,69.0589371,69.0515366,69.0469208,69.0541611,69.0304947,69.0309525,69.0141525,69.0036545,68.9837875,68.9451981,68.9588165,68.9544449,68.9217987,69.0113754,69.0094757,68.9583817]


;cluster?

;rawant=[264.6865845,264.6813354,264.6925354,264.6886597,264.6877747,264.6783752,264.6772461,264.6690979]

;decwant=[69.0454559,69.0438309,69.0457993,69.0440140,69.0423203,69.0475464,69.0443039,69.0433273]

;rawant=[265.28787]
;decwant=[69.04287]

;print, rawant, decwant

best = 0; do find only the best match, 0 = find all nearest matches

;rawant=[ 265.2015027, 264.9978887, 264.8761347, 265.1247545, 264.9859035, 265.0792373, 265.0322105, 265.0629030, 265.0369870, 265.0909960, 265.0188216, 265.0274689, 264.9765675, 264.9657073, 265.0337382, 265.0078437, 265.0754356, 265.0513790, 264.9954781, 264.9817260, 265.1850693, 265.0873221, 265.2482825, 264.9339984, 264.9297584, 264.8658059, 264.9658749, 265.0519856]

; decwant=[68.9390454, 68.9462142, 68.9551668, 68.9603397, 68.9651135, 68.9739565, 68.9816376, 68.9827230, 68.9855146, 68.9867463, 68.9870191, 68.9881279, 68.9883540, 68.9946784, 68.9974074, 68.9985320, 69.0030209, 69.0040907, 69.0056003, 69.0056641, 69.0112103, 69.0135094, 69.0186755, 69.0195889, 69.0316157, 69.0314458, 69.0333570, 69.0507640]

;openw, outlun, '/Users/jkrick/nep/pop3/want.reg', /get_lun
;printf, outlun, 'fk5'
;for w=0, n_elements(rawant)-1 do printf, outlun, 'circle( ', rawant(w), decwant(w), ' 3")'
;close, outlun
;free_lun, outlun 

;openw, outlun, '/Users/jkrick/nep/have.reg', /get_lun
;printf, outlun, 'fk5'



;readcol, '/Users/jkrick/nep/pop3/mytransients.txt',rawant, decwant, format="A"



restore, '/Users/jkrick/idlbin/idlbin_nutella/objectnew.sav'

openw, outlun, '/Users/jkrick/Desktop/foreric_mips24.txt',/get_lun
printf, outlun, "num, ra, dec, phot z, umag, umag err, gmag, gmagerr, rmag, rmagerr, imag, imagerr, acsmag, acsmagerr, zmag, zmagerr, flamingo J mag, J err, wirc J mag, J err, wirc H mag, H err, wirc K mag, K err, irac1, irac2, irac3, irac4, mips24, mips24 err, mips 70, mips70 err"


; create initial arrays
m=n_elements(rawant)
ir=n_elements(objectnew.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching  to objectnew'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
;   print, q, rawant[q]
   print, 'working on new obj'
   dist=sphdist( rawant[q], decwant[q], objectnew.ra, objectnew.dec, /degrees )
   print, 'n_elements(dist)', n_elements(dist)
   sep=min(dist,ind)
;   print, sep
   if best eq 0 then begin
      for z = 0l, n_elements(dist) - 1 do begin
         if dist(z) lt 0.0442 then begin; 0.0008 then begin
            mmatch[q]=ind

        ;put all near-IR into AB mags
         if objectnew[z].flamjmag gt 0 and objectnew[z].flamjmag ne 99 then begin
            fab = 1594.*10^(objectnew[z].flamjmag/(-2.5))
            jmagab = -2.5*alog10(fab) +8.926
         endif else begin
            jmagab = objectnew[z].flamjmag
         endelse
         
         if objectnew[z].wircjmag gt 0 and objectnew[z].wircjmag ne 99 then begin
            wircjab = 1594.*10^(objectnew[z].wircjmag/(-2.5))
            wircjmagab = -2.5*alog10(wircjab) +8.926
         endif else begin
            wircjmagab = objectnew[z].wircjmag
         endelse
         
         if objectnew[z].wirchmag gt 0 and objectnew[z].wirchmag ne 99 then begin
            wirchab = 1024.*10^(objectnew[z].wirchmag/(-2.5))
            wirchmagab = -2.5*alog10(wirchab) +8.926
         endif else begin
            wirchmagab = objectnew[z].wirchmag
         endelse
         
         if objectnew[z].wirckmag gt 0 and objectnew[z].wirckmag ne 99 then begin
            wirckab = 666.8*10^(objectnew[z].wirckmag/(-2.5))
            wirckmagab = -2.5*alog10(wirckab) +8.926
         endif else begin
            wirckmagab = objectnew[z].wirckmag
         endelse
         

         if objectnew[z].mips24flux gt 0 then begin
            printf,outlun, z, objectnew[z].ra, objectnew[z].dec, objectnew[z].zphot, objectnew[z].umag, objectnew[z].umagerr, objectnew[z].gmag, objectnew[z].gmagerr,objectnew[z].rmag, objectnew[z].rmagerr,objectnew[z].imag, objectnew[z].imagerr, objectnew[z].acsmag, objectnew[z].acsmagerr, objectnew[z].zmagbest, objectnew[z].zmagerrbest, jmagab, objectnew[z].flamjmagerr, wircjmagab, objectnew[z].wircjmagerr, wirchmagab, objectnew[z].wirchmagerr, wirckmagab, objectnew[z].wirckmagerr, objectnew[z].irac1mag,objectnew[z].irac2mag,objectnew[z].irac3mag,objectnew[z].irac4mag, objectnew[z].mips24flux, objectnew[z].mips24fluxerr,objectnew[z].mips70flux, objectnew[z].mips70fluxerr, format = '(I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'
         endif

;

;            printf,outlun, dist(z),  z, objectnew[z].ra, objectnew[z].dec, objectnew[z].zphot, objectnew[z].umag, objectnew[z].umagerr, objectnew[z].gmag, objectnew[z].gmagerr,objectnew[z].rmag, objectnew[z].rmagerr,objectnew[z].imag, objectnew[z].imagerr, objectnew[z].acsmag, objectnew[z].acsmagerr, objectnew[z].zmagbest, objectnew[z].zmagerrbest, objectnew[z].flamjmag, objectnew[z].flamjmagerr, objectnew[z].wircjmag, objectnew[z].wircjmagerr, objectnew[z].wirchmag, objectnew[z].wirchmagerr, objectnew[z].wirckmag, objectnew[z].wirckmagerr, objectnew[z].irac1mag,objectnew[z].irac2mag,objectnew[z].irac3mag,objectnew[z].irac4mag, objectnew[z].mips24flux, objectnew[z].mips24fluxerr,objectnew[z].mips70flux, objectnew[z].mips70fluxerr, format = '(F10.5,I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'
;
         endif 
      endfor
   endif
   if best eq 1 then begin
      if (sep LT 0.001)  then begin
         ;put all near-IR into AB mags
         if objectnew[ind].flamjmag gt 0 and objectnew[ind].flamjmag ne 99 then begin
            fab = 1594.*10^(objectnew[ind].flamjmag/(-2.5))
            jmagab = -2.5*alog10(fab) +8.926
         endif else begin
            jmagab = objectnew[ind].flamjmag
         endelse
         
         if objectnew[ind].wircjmag gt 0 and objectnew[ind].wircjmag ne 99 then begin
            wircjab = 1594.*10^(objectnew[ind].wircjmag/(-2.5))
            wircjmagab = -2.5*alog10(wircjab) +8.926
         endif else begin
            wircjmagab = objectnew[ind].wircjmag
         endelse
         
         if objectnew[ind].wirchmag gt 0 and objectnew[ind].wirchmag ne 99 then begin
            wirchab = 1024.*10^(objectnew[ind].wirchmag/(-2.5))
            wirchmagab = -2.5*alog10(wirchab) +8.926
         endif else begin
            wirchmagab = objectnew[ind].wirchmag
         endelse
         
         if objectnew[ind].wirckmag gt 0 and objectnew[ind].wirckmag ne 99 then begin
            wirckab = 666.8*10^(objectnew[ind].wirckmag/(-2.5))
            wirckmagab = -2.5*alog10(wirckab) +8.926
         endif else begin
            wirckmagab = objectnew[ind].wirckmag
         endelse
         

            print,ind, objectnew[ind].ra, objectnew[ind].dec, objectnew[ind].zphot, objectnew[ind].umag, objectnew[ind].umagerr, objectnew[ind].gmag, objectnew[ind].gmagerr,objectnew[ind].rmag, objectnew[ind].rmagerr,objectnew[ind].imag, objectnew[ind].imagerr, objectnew[ind].acsmag, objectnew[ind].acsmagerr, objectnew[ind].zmagbest, objectnew[ind].zmagerrbest, jmagab, objectnew[ind].flamjmagerr, wircjmagab, objectnew[ind].wircjmagerr, wirchmagab, objectnew[ind].wirchmagerr, wirckmagab, objectnew[ind].wirckmagerr, objectnew[ind].irac1mag,objectnew[ind].irac2mag,objectnew[ind].irac3mag,objectnew[ind].irac4mag, objectnew[ind].mips24flux, objectnew[ind].mips24fluxerr,objectnew[ind].mips70flux, objectnew[ind].mips70fluxerr, format = '(I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'
;

         mmatch[q]=ind
;         printf, outlun, 'circle( ', objectnew[ind].ra, objectnew[ind].dec, ' 5")'
      endif
   endif

endfor

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

print,"Finished at "+systime()
print, 'enter /Users/jkrick/nutella/zphot/zphot.target.param'
;plothyperz, mmatch, '/Users/jkrick/nutella/nep/junk.ps'
;/Users/jkrick/zphot/zphot.target.param

close, outlun
free_lun, outlun 

return, mmatch
END

;for j= 0, n_elements(rawant) - 1 do begin
;   print, rawant(j)
;   for i =0l, n_elements(objectnew.rmaga) - 1 do begin

 ;    if objectnew[i].ra gt rawant(j)-deltara and  objectnew[i].ra lt rawant(j) + deltara then begin
;         if objectnew[i].dec GT decwant(j)-deltadec and  objectnew[i].dec lt decwant(j)+deltadec then begin
 
;      if objectnew[i].zphot eq 0.9449 and objectnew[i].mass gt 1E11 then begin
;         keeper(j) = i
 
;          print, i, objectnew[i].ra, objectnew[i].dec, objectnew[i].zphot, objectnew[i].umaga, objectnew[i].umagerra, objectnew[i].gmaga, objectnew[i].gmagerra,objectnew[i].rmaga, objectnew[i].rmagerra,objectnew[i].imaga, objectnew[i].imagerra, objectnew[i].flamjmag, objectnew[i].flamjmagerr, objectnew[i].wircjmag, objectnew[i].wircjmagerr, objectnew[i].wirchmag, objectnew[i].wirchmagerr, objectnew[i].wirckmag, objectnew[i].wirckmagerr, objectnew[i].irac1mag,objectnew[i].irac2mag,objectnew[i].irac3mag,objectnew[i].irac4mag, objectnew[i].mips24flux, objectnew[i].mips24fluxerr,format = '(I10,F10.5,F10.6,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'

  
;endfor

;variable = objectnew[keeper]
;print, variable[3].rmaga
;save, variable, filename='/Users/jkrick/idlbin/variable.sav'

;ashby
;rawant = [264.85087500,264.87658333 ,265.10229167 ,265.12950000 ]    
;decwant = [69.05496944,69.03311944,68.88478889 ,69.05895278 ]

;jason
;rawant=[265.0628052,265.0905457,265.0192566,265.0273743,264.9657898,265.0337219,265.0077820,265.0514832,264.9815063]

; decwant=[68.9827576, 68.9866257, 68.9869919, 68.9880600, 68.9945679, 68.9974365, 68.9985123, 69.0042648, 69.0056915]

;zwant = 0.94
;keeper = fltarr(n_elements(rawant))
;deltara= 3.0E-4
;deltadec=3.0E-4
