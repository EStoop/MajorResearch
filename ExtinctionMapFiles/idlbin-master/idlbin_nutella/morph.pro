pro morph

restore, '/Users/jkrick/idlbin/object.sav'

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

;---------------------------------------------------------

ps_open, filename='/Users/jkrick/hst/morph2.ps',/portrait,/square,/color

;what are the extreme outliers?
;seem to be CR on edge of frame, or
;are stars with CR at their centers which did not get rejected?

for i = 0l, n_elements(object.ra) - 1, 1 do begin
;   if object[i].acsgini lt 0 then begin
;      print, "g lt 0"
;      print, x(i),y(i),mag(i),s(i)
;   endif
;   if c(i) lt 0 then begin
;      print, "c lt 0"
;      print, x(i),y(i),mag(i),s(i)
;   endif

;don't include stars, or questionable measurements

;   if object[i].s(i) gt 0.8 then c(i) =-100
   if object[i].acsmag gt 27.0 then object[i].acsconcentration = -100
   if object[i].acsellip gt 0.9 then object[i].acsconcentration = -100
;   if object[i]flag(i) gt 0 then c(i) =  -100
   ;only bright galaxies?
;   if mag(i) gt 27.0 then c(i)=-100
   if object[i].acsconcentration le 0 then object[i].acsconcentration = -100
   if object[i].acsgini le 0 then object[i].acsgini = -100

;aim for a phot-z range
;   if object[i].zphot gt 0.3 then object[i].acsconcentration=-100

endfor


;color code by surface brightness

plot, object.acsconcentration , object.acsgini, psym = 3, xthick=3,ythick=3,charthick=3, xtitle='concentration', ytitle='gini coefficient',xrange=[0,1], yrange=[0,1]

sb1 = where(object.acscentralmu lt 29 and object.acscentralmu gt 15 )
sb2 = where(object.acscentralmu lt 30 and object.acscentralmu ge 29 )
sb3 = where(object.acscentralmu lt 30.5 and object.acscentralmu ge 30 )
sb4 = where(object.acscentralmu ge 30.5 )

oplot, object[sb1].acsconcentration, object[sb1].acsgini, psym = 2, color = redcolor
oplot, object[sb2].acsconcentration, object[sb2].acsgini, psym = 2, color = greencolor
oplot, object[sb3].acsconcentration, object[sb3].acsgini, psym = 2, color = bluecolor
oplot, object[sb4].acsconcentration, object[sb4].acsgini, psym = 2, color = purplecolor


;-------------------------------------------------------------------
;fitting section
;---------------------------------------------------------------------

;sort them in concentration space to fit by bins
sortindex= sort(object.acsconcentration)
csort= object[sortindex].acsconcentration
gsort = object[sortindex].acsgini


c1 = fltarr(n_elements(object.ra))
g1 = fltarr(n_elements(object.ra))
c2 = fltarr(n_elements(object.ra))
g2 = fltarr(n_elements(object.ra))
t= 0
u = 0
;make the cutoff at c = 0.2, fit before 0.2 and between 0.2 and 0.6  with two
;different functions

FOR se = 0l, n_elements(object.ra)-1, 1 DO begin
   IF (object[se].acsconcentration LE 0.2) AND (object[se].acsconcentration GT 0.0)  THEN BEGIN
      c1(t) = object[se].acsconcentration
      g1(t) = object[se].acsgini
      t = t +1
   ENDIF 
   IF (object[se].acsconcentration LT 0.6) AND (object[se].acsconcentration GT 0.2) THEN BEGIN
      c2(u) = object[se].acsconcentration
      g2(u) = object[se].acsgini
      u = u + 1
   ENDIF

ENDFOR

print, t, " the number of concentrations less than 0.2"
c1 = c1[0:t-1]
g1 = g1[0:t-1]
c2 = c2[0:u-1]
g2 = g2[0:u-1]

;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( mgalaxy.rmag, mgalaxy.vr, yfit, sig, coeff_sig)
coeff1 = ROBUST_LINEFIT(c1,g1, yfit1, sig1, coeff_sig1)
coeff2 = ROBUST_LINEFIT(c2,g2, yfit2, sig2, coeff_sig2)
print, "sig1, sig2", sig1,sig2

err = dindgen(n_elements(c1) - 1) - dindgen(n_elements(c1) - 1) + 1
start = [1.0,0.1]
;c3=[c1,c2]
;yfit3 = [yfit1,yfit2]

sortindex = sort(c2)
sortc2 = c2[sortindex]
sortyfit2 = yfit2[sortindex]
result = MPFITFUN('linear',sortc2, sortyfit2,err, start)
print, "result", result
xvals = findgen(11) /10
oplot, xvals, ((result(0))*xvals) + result(1), thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)+ sig1 , thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)- sig1 , thick = 3;, color = colors.orange


;color code by magnitude

plot, object.acsconcentration , object.acsgini, psym = 3, xthick=3,ythick=3,charthick=3, xtitle='concentration', ytitle='gini coefficient',xrange=[0,1], yrange=[0,1]

for count =0l,n_elements(object.ra) - 1, 1 do begin
   if object[count].acsmag lt 26 and object[count].acsmag gt 15 then begin
      xyouts, object[count].acsconcentration, object[count].acsgini, 'x', color = redcolor,alignment=0.5
   endif
   if object[count].acsmag lt 27 and object[count].acsmag ge 26 then begin
      xyouts, object[count].acsconcentration, object[count].acsgini, 'x', color = greencolor,alignment=0.5
   endif
   if object[count].acsmag lt 28 and object[count].acsmag ge 27 then begin
      xyouts, object[count].acsconcentration, object[count].acsgini, 'x', color = bluecolor,alignment=0.5
   endif
   if  object[count].acsmag ge 28 then begin
      xyouts, object[count].acsconcentration, object[count].acsgini , 'x',alignment=0.5
   endif

endfor



ps_close, /noprint,/noid

undefine, object
end

;;----------------------------------------------------------
;;put the morpheus values into the catalog.

;ftab_ext, '/Users/jkrick/hst/raw/wholeacs_catalog.fits', 'Concentration,Gini,mag_auto,SurfaceBrightness,x_image,y_image,class_star,a_image,b_image', c,g,mag,mu,x,y,s,a,b
;ftab_ext, '/Users/jkrick/hst/raw/wholeacs_catalog.fits', 'flags,x_world,y_world,isoarea_image,theta_world,signaltonoiseratio,fluxcounts,sky,skyrms',flag,ra,dec,isoarea,theta,snr,flux, sky,skyrms
;ftab_ext, '/Users/jkrick/hst/raw/wholeacs_catalog.fits', 'centralsurfacebrightness,rotationalasymmetry,halflightradiusinpixels,m20,aperturemagnitude,isophotalmagnitude,magerr_auto',cmu,asym,hlr,m20,apmag,isomag,magerrauto

;can run another version of sextractacscatalog  which will output things like fwhm or whatever else
;but won't need to be run through analyzecatalog  so the output being wrong won't matter

;autoarr = dblarr(n_elements(c))
;bestarr = autoarr
;aperarr = autoarr
;isoarr = autoarr
;newindex = 0l

; create initial arrays
;m=n_elements(c)
;ir=n_elements(object.ra)

;irmatch=fltarr(ir)
;mmatch=fltarr(m)
;irmatch[*]=-999
;mmatch[*]=-999

;print, m
;print,'Matching morpheus to object'
;print,"Starting at "+systime()
;dist=irmatch
;dist[*]=0

;for q=0l,m-1 do begin

;   dist=sphdist( ra(q), dec(q),object.acsra,object.acsdec,/degrees)
;   sep=min(dist,ind)
;   if (sep LT 0.0008)  then begin
;      mmatch[q]=ind
;;      object[ind].acsra=ra(q)
;;      object[ind].acsdec=dec(q)
;;      object[ind].acsisoarea=isoarea(q)
;;      object[ind].acsellip=1-b(q)/a(q)
;;      object[ind].acstheta=theta(q)
;      object[ind].acssnr=snr(q)
;      object[ind].acsmag=mag(q)
;      object[ind].acsflux=flux(q)
;      object[ind].acssky=sky(q)
;      object[ind].acsskyrms=skyrms(q)
;      object[ind].acsmu=mu(q)
;      object[ind].acscentralmu=cmu(q)
;      object[ind].acsgini=g(q)
;      object[ind].acsassym=asym(q)
;      object[ind].acshlightrad=hlr(q)
;      object[ind].acsconcentration=c(q)
;      object[ind].acsm20=m20(q);
;      object[ind].acsflag = flag(q)
;      object[ind].acsmagerr = magerrauto(q)

;      autoarr(newindex) = mag(q)
;      bestarr(newindex) = object[ind].acsmag
;      aperarr(newindex) = apmag(q)
;      isoarr(newindex) = isomag(q)
;      newindex = newindex + 1
;   endif 
;endfor

;print,"Finished at "+systime()
;save, object, filename='/Users/jkrick/idlbin/object.sav'

;plothist, mag, xhist, yhist, /noprint,/nan,bin=0.05,xrange=[25,30]
;plothist,cmu,xhist,yhist,/noprint,/nan,bin=0.05,xrange=[25,30]
;plothist,object.acscentralmu,xhist,yhist,/noprint,/nan,bin=0.05,xrange=[25,30]
;plothist, object.zphot, xhist, yhist, /noprint

;plot, autoarr, autoarr - bestarr, psym = 3, thick = 3, xtitle="mag_auto", ytitle = "mag_auto - mag_best", xrange=[12,30], yrange=[-6,2]
;plot, isoarr, isoarr - bestarr, psym = 3, thick = 3, xtitle="mag_iso", ytitle = "mag_iso - mag_best", xrange=[12,30], yrange=[-6,2]

;for count =0l,n_elements(object.ra) - 1, 1 do begin
;   if object[count].acscentralmu lt 29 and object[count].acscentralmu gt 15 then begin
;      xyouts, object[count].acsconcentration, object[count].acsgini, 'x', color = redcolor,alignment=0.5
;   endif
;   if object[count].acscentralmu lt 30 and object[count].acscentralmu ge 29 then begin
;      xyouts, object[count].acsconcentration, object[count].acsgini, 'x', color = greencolor,alignment=0.5
;   endif
;   if object[count].acscentralmu lt 30.5 and object[count].acscentralmu ge 30 then begin
;      xyouts, object[count].acsconcentration, object[count].acsgini, 'x', color = bluecolor,alignment=0.5
;   endif
;   if  object[count].acscentralmu ge 30.5 then begin
;      xyouts, object[count].acsconcentration, object[count].acsgini, 'x',alignment=0.5
;   endif

;endfor

