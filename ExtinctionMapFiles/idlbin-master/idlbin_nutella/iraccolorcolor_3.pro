pro iraccolorcolor
close,/all

restore, '/Users/jkrick/idlbin/objectnew.sav'

ps_open, filename='/Users/jkrick/nep/colorcolor_old.ps',/portrait,/square,/color
;redcolor = FSC_COLOR("Red", !D.Table_Size-2)
;bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
;greencolor = FSC_COLOR("Green", !D.Table_Size-4)
;yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
;cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
;orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
;purplecolor = FSC_COLOR("purple", !D.Table_Size-8)


!p.thick=3
!x.thick=3
!y.thick=3
!p.charsize=1.5
!p.charthick=3
!p.multi = [0, 0, 1]
 
vsym, /polygon, /fill

 ix=where(objectnew.irac1fluxold gt 0. and objectnew.irac4fluxold gt 0.0 and objectnew.irac3fluxold gt 00.0 and objectnew.irac2fluxold gt 0.0   and objectnew.chisq lt 50)
; ix=where(objectnew.irac1fluxold gt 0. and objectnew.irac4fluxold gt 0.0 and objectnew.irac3fluxold gt 00.0 and objectnew.irac2fluxold gt 0.0  )
print, "n_elenents(ix)", n_elements(ix)

agn = where(objectnew[ix].spt eq 6 or objectnew[ix].spt eq 7 or objectnew[ix].spt eq 13$
or objectnew[ix].spt eq 14 or objectnew[ix].spt eq 15)
print, "n_elememts(agn)", n_elements(agn)

star = where(objectnew[ix].acsclassstar eq 1 or objectnew[ix].uclassstar eq 1 or objectnew[ix].gclassstar eq 1 or objectnew[ix].rclassstar eq 1 or objectnew[ix].iclassstar eq 1 )

print, "n_elements(star)", n_elements(star)

composite = where(objectnew[ix].spt eq 5 )

starburst = where(objectnew[ix].spt eq 1 or objectnew[ix].spt eq 4)

sp = where(objectnew[ix].spt eq 10 or objectnew[ix].spt eq 11 or objectnew[ix].spt eq 12 or objectnew[ix].spt eq 9) 

ell = where(objectnew[ix].spt eq 2 or objectnew[ix].spt eq  3 );or objectnew[ix].spt eq 8)

sf = [ sp]

print, n_elements(sf), "n_elements(sf)"
print, n_elements(composite), "n_elements(composite)"

;---------------------
plot, alog10(objectnew[ix].irac3fluxold / objectnew[ix].irac1fluxold), alog10(objectnew[ix].irac4fluxold / objectnew[ix].irac2fluxold), psym = 8, xrange = [-1.2, 1.], yrange=[-1.0, 1.0],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', symsize=0.3

oplot, alog10(objectnew[ix[ell]].irac3fluxold / objectnew[ix[ell]].irac1fluxold), alog10(objectnew[ix[ell]].irac4fluxold / objectnew[ix[ell]].irac2fluxold), psym = 5, symsize=0.5, color = yellowcolor

oplot, alog10(objectnew[ix[sf]].irac3fluxold / objectnew[ix[sf]].irac1fluxold), alog10(objectnew[ix[sf]].irac4fluxold / objectnew[ix[sf]].irac2fluxold), psym = 5, symsize=0.5, color = bluecolor

oplot, alog10(objectnew[ix[agn]].irac3fluxold / objectnew[ix[agn]].irac1fluxold), alog10(objectnew[ix[agn]].irac4fluxold / objectnew[ix[agn]].irac2fluxold), psym = 5, symsize=0.3, color = redcolor


oplot, alog10(objectnew[ix[star]].irac3fluxold / objectnew[ix[star]].irac1fluxold), alog10(objectnew[ix[star]].irac4fluxold / objectnew[ix[star]].irac2fluxold), psym = 5, symsize=0.3, color = greencolor

;oplot, alog10(objectnew[ix].irac3fluxold / objectnew[ix].irac1fluxold), alog10(objectnew[ix].irac4fluxold / objectnew[ix].irac2fluxold), psym = 8, symsize=0.3


;oplot, alog10(objectnew[ix[composite]].irac3fluxold / objectnew[ix[composite]].irac1fluxold), alog10(objectnew[ix[composite]].irac4fluxold / objectnew[ix[composite]].irac2fluxold), psym = 5, symsize=0.3, color = purplecolor


oplot, [-0.2, 1.0], [-0.2,-0.2], linestyle = 2, thick =5
oplot, [-0.2, -0.2], [-0.2,0.3], linestyle = 2, thick =5
oplot, [-0.2, 1.0], [0.3,1.3], linestyle = 2,thick = 5

;--------------------
plot, objectnew[ix].irac3magold  - objectnew[ix].irac4magold, objectnew[ix].irac1magold - objectnew[ix].irac2magold, xrange = [-1.5, 2.0], yrange=[-0.8,0.7], psym = 8, /xst, /yst, xtitle='[5.8] - [8.0] (AB)', ytitle = '[3.6] - [4.5] (AB)', symsize=0.3

oplot, objectnew[ix[ell]].irac3magold  - objectnew[ix[ell]].irac4magold, objectnew[ix[ell]].irac1magold - objectnew[ix[ell]].irac2magold, psym = 5, symsize=0.3, color=yellowcolor

oplot, objectnew[ix[sf]].irac3magold  - objectnew[ix[sf]].irac4magold, objectnew[ix[sf]].irac1magold - objectnew[ix[sf]].irac2magold, psym = 5, symsize=0.3, color=bluecolor

oplot, objectnew[ix[agn]].irac3magold  - objectnew[ix[agn]].irac4magold, objectnew[ix[agn]].irac1magold - objectnew[ix[agn]].irac2magold, psym = 5, symsize=0.3, color=redcolor

oplot, objectnew[ix[star]].irac3magold  - objectnew[ix[star]].irac4magold, objectnew[ix[star]].irac1magold - objectnew[ix[star]].irac2magold, psym = 5, symsize=0.3, color=greencolor

oplot, [-0.1, -0.1], [-0.2,1.0], linestyle = 2, thick =5
oplot, [-0.1, 0.9], [-0.2,0], linestyle = 2, thick =5
oplot, [0.9, 1.3], [0.0,1.0], linestyle = 2,thick = 5

;-------------------------------------------------------
;look at the variance in the x axis of the lacy plot
mmm,  alog10(objectnew[ix].irac3fluxold / objectnew[ix].irac1fluxold), modeold, sigmaold
varold = sigmaold^2
;plothist, alog10(objectnew[ix].irac3fluxold / objectnew[ix].irac1fluxold)
print, 'mode ,sigma old', modeold, sigmaold

mmm,  alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), mode, sigma
var = sigma^2
;plothist, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux)
print, 'mode ,sigma', mode, sigma


if var gt varold then f = var / varold
if var lt varold then f = varold / var

prob = mpftest(f, n_elements(ix)-1, n_elements(ix)-1)

print, 'f, prob lacy x', f, prob
;-------------------------------------------------------
;look at the variance in the y axis of the lacy plot
mmm,  alog10(objectnew[ix].irac4fluxold / objectnew[ix].irac2fluxold), modeold, sigmaold
varold = sigmaold^2
;plothist, alog10(objectnew[ix].irac3fluxold / objectnew[ix].irac1fluxold)
print, 'mode ,sigma old', modeold, sigmaold

mmm,  alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux), mode, sigma
var = sigma^2
;plothist, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux)
print, 'mode ,sigma', mode, sigma


if var gt varold then f = var / varold
if var lt varold then f = varold / var

prob = mpftest(f, n_elements(ix)-1, n_elements(ix)-1)

print, 'f, prob lacy y', f, prob
;-------------------------------------------------------
;look at the variance in the x axis of the stern plot
mmm,   objectnew[ix].irac3magold  - objectnew[ix].irac4magold, modeold, sigmaold
varold = sigmaold^2
print, 'mode ,sigma old', modeold, sigmaold

mmm,   objectnew[ix].irac3mag  - objectnew[ix].irac4mag, mode, sigma
var = sigma^2
print, 'mode ,sigma', mode, sigma


if var gt varold then f = var / varold
if var lt varold then f = varold / var

prob = mpftest(f, n_elements(ix)-1, n_elements(ix)-1)

print, 'f, prob stern x', f, prob

;-------------------------------------------------------
;look at the variance in the y axis of the stern plot
mmm,   objectnew[ix].irac1magold  - objectnew[ix].irac2magold, modeold, sigmaold
varold = sigmaold^2
print, 'mode ,sigma old', modeold, sigmaold

mmm,   objectnew[ix].irac1mag  - objectnew[ix].irac2mag, mode, sigma
var = sigma^2
print, 'mode ,sigma', mode, sigma


if var gt varold then f = var / varold
if var lt varold then f = varold / var

prob = mpftest(f, n_elements(ix)-1, n_elements(ix)-1)

print, 'f, prob stern x', f, prob




ps_close, /noprint,/noid



;---------------------
;test
;what is the redshift of the AGn to the left of the Lacy wedge?

test =  where(alog10(objectnew[ix[agn]].irac3fluxold / objectnew[ix[agn]].irac1fluxold) lt -0.2)
print, 'test', n_elements(test)

plothist, objectnew[ix[agn[test]]].zphot, bin = 0.1

end

