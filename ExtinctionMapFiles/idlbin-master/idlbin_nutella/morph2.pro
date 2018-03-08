pro morph

restore, '/Users/jkrick/idlbin/object.sav'



ps_open, filename='/Users/jkrick/hst/morph2.ps',/portrait,/square,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

vsym, /polygon, /fill
!p.symsize=0.3
;---------------------------------------------------------------------
;color code by surface brightness
;---------------------------------------------------------------------

plot, object.acsconcentration , object.acsgini, psym = 3, xthick=3,ythick=3,charthick=3, xtitle='concentration', ytitle='gini coefficient',xrange=[0,1], yrange=[0,1]

sb1 = where(object.acsmag lt 27. and object.acsclassstar lt 1. and object.acscentralmu lt 29 and object.acscentralmu gt 15 )
sb2 = where(object.acsmag lt 27. and object.acsclassstar lt 1. and object.acscentralmu lt 30 and object.acscentralmu ge 29 )
sb3 = where(object.acsmag lt 27. and object.acsclassstar lt 1. and object.acscentralmu lt 30.5 and object.acscentralmu ge 30 )
sb4 = where(object.acsmag lt 27. and object.acsclassstar lt 1. and object.acscentralmu ge 30.5 )

oplot, object[sb1].acsconcentration, object[sb1].acsgini, psym = 8, color = redcolor, symsize=0.3
oplot, object[sb2].acsconcentration, object[sb2].acsgini, psym = 8, color = greencolor, symsize=0.3
oplot, object[sb3].acsconcentration, object[sb3].acsgini, psym = 8, color = bluecolor, symsize=0.3
oplot, object[sb4].acsconcentration, object[sb4].acsgini, psym = 8, color = purplecolor, symsize=0.3

;-------------------------------------------------------------------
;fitting section
;---------------------------------------------------------------------

;sort them in concentration space to fit by bins
sortindex= sort(object.acsconcentration)
csort= object[sortindex].acsconcentration
gsort = object[sortindex].acsgini

;make the cutoff at c = 0.2, fit before 0.2 and between 0.2 and 0.6  with two
;different functions

c1 = where(object.acsconcentration le 0.2 and object.acsconcentration gt 0)
c2 = where(object.acsconcentration le 0.6 and object.acsconcentration gt 0.2)



print, n_elements(c1), " the number of concentrations less than 0.2"

;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( mgalaxy.rmag, mgalaxy.vr, yfit, sig, coeff_sig)
coeff1 = ROBUST_LINEFIT(object[c1].acsconcentration,object[c1].acsgini, yfit1, sig1, coeff_sig1)
coeff2 = ROBUST_LINEFIT(object[c2].acsconcentration,object[c2].acsgini, yfit2, sig2, coeff_sig2)
print, "sig1, sig2", sig1,sig2

err = dindgen(n_elements(c1) - 1) - dindgen(n_elements(c1) - 1) + 1
start = [1.0,0.1]
;c3=[c1,c2]
;yfit3 = [yfit1,yfit2]

sortindex = sort(object[c2].acsconcentration)
sortc2 = object[c2[sortindex]].acsconcentration
sortyfit2 = yfit2[sortindex]
result = MPFITFUN('linear',sortc2, sortyfit2,err, start)
print, "result", result
xvals = findgen(11) /10
oplot, xvals, ((result(0))*xvals) + result(1), thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)+ sig1 , thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)- sig1 , thick = 3;, color = colors.orange

;---------------------------------------------------------------------
;color code by magnitude
;---------------------------------------------------------------------

plot, object.acsconcentration , object.acsgini, psym = 3, xthick=3,ythick=3,charthick=3, xtitle='concentration', ytitle='gini coefficient',xrange=[0,1], yrange=[0,1]

mag1 = where(object.acsmag lt 26 and object.acsmag gt 15 and object.acsclassstar lt 1.0)
mag2 = where(object.acsmag lt 27 and object.acsmag ge 26 and object.acsclassstar lt 1.0)
mag3 = where(object.acsmag lt 28 and object.acsmag ge 27 and object.acsclassstar lt 1.0)
mag4 = where(object.acsmag ge 28)

oplot, object[mag1].acsconcentration, object[mag1].acsgini, psym = 8, color = redcolor, symsize=0.3
oplot, object[mag2].acsconcentration, object[mag2].acsgini, psym = 8, color = greencolor, symsize=0.3
oplot, object[mag3].acsconcentration, object[mag3].acsgini, psym = 8, color = bluecolor, symsize=0.3
oplot, object[mag4].acsconcentration, object[mag4].acsgini, psym = 8, color = purplecolor, symsize=0.3


;ell = where(object.spt eq 1 or object.spt eq  2 or object.spt eq 8 and object.prob gt 0.9 )
;sp = where(object.spt eq 9 or object.spt eq 10 or object.spt eq 11 or object.spt eq 12 and object.prob gt 0.9) 
;xyouts, 0.7, 0.04, 'hyperz spirals', color = greencolor
;xyouts, 0.7, 0.07, 'hyperz ellipticals', color = redcolor

;---------------------------------------------------------------------
;compare to COSMOS swarms, ie Scarlata et al.
;---------------------------------------------------------------------


good = where(object.acsmag lt 24. and object.acsmag gt 0  and object.acsclassstar lt 1.0 and object.irac1mag gt 0 and object.irac1mag lt 90 )
print, n_elements(good), "   n_elements(good)", stddev(object[good].acsm20)

;M20 vs gini
plot, alog10(object[good].acsm20) , object[good].acsgini, psym = 8, xthick=3,ythick=3,charthick=3, xtitle='m20', ytitle='gini coefficient',xrange=[0,-3], yrange=[0,1]
im_hessplot, alog10(object[good].acsm20) , object[good].acsgini,  xtitle='m20', ytitle='gini coefficient',xrange=[0,-3], yrange=[0,1]

;plot M20 vs C
plot, alog10(object[good].acsm20), object[good].acsconcentration, psym = 8, xthick=3,ythick=3,charthick=3, xtitle='M20', ytitle='concentration', xrange=[1,-3],yrange=[0,1]
im_hessplot, alog10(object[good].acsm20), object[good].acsconcentration, xtitle='M20', ytitle='concentration', xrange=[0,-3],yrange=[0,1]

;plot A vs C
plot, object[good].acsassym, object[good].acsconcentration, psym = 8, xthick=3,ythick=3,charthick=3, xtitle='assymetry', ytitle='concentration',xrange=[0,0.6], yrange=[0,1]
im_hessplot, object[good].acsassym, object[good].acsconcentration,  xtitle='assymetry', ytitle='concentration',xrange=[0,0.6], yrange=[0,1]


ps_close, /noprint,/noid


end
