pro iraccolorcolor
close,/all

restore, '/Users/jkrick/idlbin/objectnew.sav'

ps_open, filename='/Users/jkrick/nep/colorcolor.ps',/portrait,/square,/color
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
plot, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux), psym = 8, xrange = [-1.2, 1.], yrange=[-1.0, 1.0],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', symsize=0.3

oplot, alog10(objectnew[ix[ell]].irac3flux / objectnew[ix[ell]].irac1flux), alog10(objectnew[ix[ell]].irac4flux / objectnew[ix[ell]].irac2flux), psym = 5, symsize=0.5, color = yellowcolor

oplot, alog10(objectnew[ix[sf]].irac3flux / objectnew[ix[sf]].irac1flux), alog10(objectnew[ix[sf]].irac4flux / objectnew[ix[sf]].irac2flux), psym = 5, symsize=0.5, color = bluecolor

oplot, alog10(objectnew[ix[agn]].irac3flux / objectnew[ix[agn]].irac1flux), alog10(objectnew[ix[agn]].irac4flux / objectnew[ix[agn]].irac2flux), psym = 5, symsize=0.3, color = redcolor


oplot, alog10(objectnew[ix[star]].irac3flux / objectnew[ix[star]].irac1flux), alog10(objectnew[ix[star]].irac4flux / objectnew[ix[star]].irac2flux), psym = 5, symsize=0.3, color = greencolor

;oplot, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux), psym = 8, symsize=0.3


;oplot, alog10(objectnew[ix[composite]].irac3flux / objectnew[ix[composite]].irac1flux), alog10(objectnew[ix[composite]].irac4flux / objectnew[ix[composite]].irac2flux), psym = 5, symsize=0.3, color = purplecolor


oplot, [-0.2, 1.0], [-0.2,-0.2], linestyle = 2, thick =5
oplot, [-0.2, -0.2], [-0.2,0.3], linestyle = 2, thick =5
oplot, [-0.2, 1.0], [0.3,1.3], linestyle = 2,thick = 5

;--------------------
plot, objectnew[ix].irac3mag  - objectnew[ix].irac4mag, objectnew[ix].irac1mag - objectnew[ix].irac2mag, xrange = [-1.5, 2.0], yrange=[-0.8,0.7], psym = 8, /xst, /yst, xtitle='[5.8] - [8.0] (AB)', ytitle = '[3.6] - [4.5] (AB)', symsize=0.3

oplot, objectnew[ix[ell]].irac3mag  - objectnew[ix[ell]].irac4mag, objectnew[ix[ell]].irac1mag - objectnew[ix[ell]].irac2mag, psym = 5, symsize=0.3, color=yellowcolor

oplot, objectnew[ix[sf]].irac3mag  - objectnew[ix[sf]].irac4mag, objectnew[ix[sf]].irac1mag - objectnew[ix[sf]].irac2mag, psym = 5, symsize=0.3, color=bluecolor

oplot, objectnew[ix[agn]].irac3mag  - objectnew[ix[agn]].irac4mag, objectnew[ix[agn]].irac1mag - objectnew[ix[agn]].irac2mag, psym = 5, symsize=0.3, color=redcolor

oplot, objectnew[ix[star]].irac3mag  - objectnew[ix[star]].irac4mag, objectnew[ix[star]].irac1mag - objectnew[ix[star]].irac2mag, psym = 5, symsize=0.3, color=greencolor

oplot, [-0.1, -0.1], [-0.2,1.0], linestyle = 2, thick =5
oplot, [-0.1, 0.9], [-0.2,0], linestyle = 2, thick =5
oplot, [0.9, 1.3], [0.0,1.0], linestyle = 2,thick = 5

ps_close, /noprint,/noid



;---------------------
;test
;what is the redshift of the AGn to the left of the Lacy wedge?

test =  where(alog10(objectnew[ix[agn]].irac3flux / objectnew[ix[agn]].irac1flux) lt -0.2)
print, 'test', n_elements(test)

;plothist, objectnew[ix[agn[test]]].zphot, bin = 0.1

end

