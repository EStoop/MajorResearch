pro iraccolorcolor
close,/all

restore, '/Users/jkrick/idlbin/objectnew.sav'

ps_open, filename='/Users/jkrick/nep/colorcolor.ps',/portrait,/square,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)


!p.thick=3
!x.thick=3
!y.thick=3
!p.charsize=1.5
!p.charthick=3
!p.multi = [0, 0, 1]
 
vsym, /polygon, /fill
;plot, objectnew.irac3mag - objectnew.irac4mag, objectnew.irac1mag - objectnew.irac2mag, psym = 8, xrange = [-1.2, 2.9], yrange=[-0.75, 1.0],/xst,/yst,symsize=0.4, xtitle='ch3 - ch4', ytitle = 'ch1 - ch2'
; oplot,[1,1]*(-0.1),[-0.2,1]
; oplot,[-0.1,1.0],[-0.2,0]
; oplot,[1.0,1.3],[0.0,1.0]
; oplot,[1,1]*(-0.1),[-0.2,1]

 ix=where(objectnew.irac1mag lt 90.0 and objectnew.irac1mag gt 0.0 and objectnew.irac4flux gt 0.0 and objectnew.irac4mag lt 90. and objectnew.irac3mag gt 0.0 and objectnew.irac3mag lt 90. and objectnew.irac2mag gt 0.0 and objectnew.irac2mag lt 90. )
print, "n_elenents(ix)", n_elements(ix)

plot, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux), psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.5],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', ticklen = 1, symsize=0.3

;im_hessplot, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux),  xrange = [-1.5, 1.], yrange=[-1.0, 1.5],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', ticklen = 1, nbin2d=60
;-----------------------------------------------------
plot, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux), psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.5],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', ticklen = 1, symsize=0.3

mips24 = where(objectnew[ix].mips24mag gt 0 and objectnew[ix].mips24mag lt 90)
oplot,  alog10(objectnew[ix[mips24]].irac3flux / objectnew[ix[mips24]].irac1flux), alog10(objectnew[ix[mips24]].irac4flux / objectnew[ix[mips24]].irac2flux), psym = 8,symsize=0.3, Color=redcolor

xyouts, -1.2, 1.2, '24 micron sources', Color  = redcolor
;-----------------------------------------------------

plot, alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux), alog10(objectnew[ix].irac4flux / objectnew[ix].irac2flux), psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.5],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', ticklen = 1, symsize=0.3

mips70 = where(objectnew[ix].mips70flux gt 0 )
oplot,  alog10(objectnew[ix[mips70]].irac3flux / objectnew[ix[mips70]].irac1flux), alog10(objectnew[ix[mips70]].irac4flux / objectnew[ix[mips70]].irac2flux), psym = 8,symsize=0.3, Color=redcolor

xyouts, -1.2, 1.2, '70 micron sources', Color  = redcolor
;-----------------------------------------------------

star = where(objectnew[ix].acsclassstar eq 1 or objectnew[ix].uclassstar eq 1 or objectnew[ix].gclassstar eq 1 or objectnew[ix].rclassstar eq 1 or objectnew[ix].iclassstar eq 1 )

agn = where(objectnew[ix].spt eq 7 or objectnew[ix].spt eq 8 or objectnew[ix].spt eq 14 or objectnew[ix].spt eq 15)
agninterest = where(objectnew[ix[agn]].prob gt 0.5 and alog10(objectnew[ix[agn]].irac3flux / objectnew[ix[agn]].irac1flux) gt 0.0 and objectnew[ix[agn]].acsmag gt 0 and objectnew[ix[agn]].acsmag lt 90)
print, n_elements(agninterest), 'agninterest'



composite = where(objectnew[ix].spt eq 4 or objectnew[ix].spt eq 6)

starburst = where(objectnew[ix].spt eq 1 or objectnew[ix].spt eq 5)
starburstinterest = where(objectnew[ix[starburst]].prob gt 0.5 and alog10(objectnew[ix[starburst]].irac3flux / objectnew[ix[starburst]].irac1flux) gt 0.0 and objectnew[ix[starburst]].acsmag gt 0 and objectnew[ix[starburst]].acsmag lt 90)
print, n_elements(starburstinterest), 'starburstinterest'

ell = where(objectnew[ix].spt eq 2 or objectnew[ix].spt eq  3 or objectnew[ix].spt eq 9)

sp = where(objectnew[ix].spt eq 10 or objectnew[ix].spt eq 11 or objectnew[ix].spt eq 12 or objectnew[ix].spt eq 13) 

test  = where(objectnew[ix].acsmag gt 0 and objectnew[ix].acsmag lt 90 and alog10(objectnew[ix].irac3flux / objectnew[ix].irac1flux) lt -1.3)
print, 'test', n_elements(test)

plot, alog10(objectnew[ix[star]].irac3flux / objectnew[ix[star]].irac1flux), alog10(objectnew[ix[star]].irac4flux / objectnew[ix[star]].irac2flux), psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.0],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', symsize=0.3 ; ticklen = 1,


oplot, alog10(objectnew[ix[ell]].irac3flux / objectnew[ix[ell]].irac1flux), alog10(objectnew[ix[ell]].irac4flux / objectnew[ix[ell]].irac2flux), psym = 8,  symsize=0.3, Color=redcolor

oplot, alog10(objectnew[ix[sp]].irac3flux / objectnew[ix[sp]].irac1flux), alog10(objectnew[ix[sp]].irac4flux / objectnew[ix[sp]].irac2flux), psym = 8,  symsize=0.3, Color=redcolor


oplot, alog10(objectnew[ix[starburst]].irac3flux / objectnew[ix[starburst]].irac1flux), alog10(objectnew[ix[starburst]].irac4flux / objectnew[ix[starburst]].irac2flux), psym = 8,  symsize=0.3, Color=redcolor

oplot, alog10(objectnew[ix[agn]].irac3flux / objectnew[ix[agn]].irac1flux), alog10(objectnew[ix[agn]].irac4flux / objectnew[ix[agn]].irac2flux), psym = 8,  symsize=0.3, Color=bluecolor

;oplot, alog10(objectnew[ix[composite]].irac3flux / objectnew[ix[composite]].irac1flux), alog10(objectnew[ix[composite]].irac4flux / objectnew[ix[composite]].irac2flux), psym = 8,  symsize=0.3, Color=yellowcolor

oplot, alog10(objectnew[ix[star]].irac3flux / objectnew[ix[star]].irac1flux), alog10(objectnew[ix[star]].irac4flux / objectnew[ix[star]].irac2flux), psym = 8,  symsize=0.3

oplot, [-0.2,1.3], [-0.2,-0.2]
oplot, [-0.2,-0.2], [-0.2,0.3]
oplot, [-0.2,1.3], [0.3,1.5]

xyouts, -1.2, 0.9, "stars"
xyouts, -1.2, 0.85, "AGN", color = bluecolor
xyouts, -1.2, 0.8, "E's and S's", color=redcolor
;xyouts, -1.2, 0.75, "starbursts", color=greencolor

;---------------------------------------------------
;try to separate in redshift space


z1e = where(objectnew[ix[ell]].zphot le 0.2 and objectnew[ix[ell]].zphot gt 0.1)
z1s = where(objectnew[ix[sp]].zphot le 0.2 and objectnew[ix[sp]].zphot gt 0.1)

z2e = where(objectnew[ix[ell]].zphot le 0.4 and objectnew[ix[ell]].zphot gt 0.2)
z2s = where(objectnew[ix[sp]].zphot le 0.4 and objectnew[ix[sp]].zphot gt 0.2)

z3e = where(objectnew[ix[ell]].zphot le 0.6 and objectnew[ix[ell]].zphot gt 0.4)
z3s = where(objectnew[ix[sp]].zphot le 0.6 and objectnew[ix[sp]].zphot gt 0.4)

z4e = where(objectnew[ix[ell]].zphot le 0.8 and objectnew[ix[ell]].zphot gt 0.6)
z4s = where(objectnew[ix[sp]].zphot le 0.8 and objectnew[ix[sp]].zphot gt 0.6)

z5e = where(objectnew[ix[ell]].zphot le 1.0 and objectnew[ix[ell]].zphot gt 0.8)
z5s = where(objectnew[ix[sp]].zphot le 1.0 and objectnew[ix[sp]].zphot gt 0.8)



plot, alog10(objectnew[ix[ell[z1e]]].irac3flux / objectnew[ix[ell[z1e]]].irac1flux), alog10(objectnew[ix[ell[z1e]]].irac4flux / objectnew[ix[ell[z1e]]].irac2flux), psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.0],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', symsize=0.3 ,/nodata, title = "ellipticals as a function of z"; ticklen = 1,

oplot, alog10(objectnew[ix[ell[z1e]]].irac3flux / objectnew[ix[ell[z1e]]].irac1flux), alog10(objectnew[ix[ell[z1e]]].irac4flux / objectnew[ix[ell[z1e]]].irac2flux), psym = 8, symsize=0.3, color=greencolor

oplot, alog10(objectnew[ix[ell[z2e]]].irac3flux / objectnew[ix[ell[z2e]]].irac1flux), alog10(objectnew[ix[ell[z2e]]].irac4flux / objectnew[ix[ell[z2e]]].irac2flux), psym = 8, symsize=0.3,color = yellowcolor
;
oplot, alog10(objectnew[ix[ell[z3e]]].irac3flux / objectnew[ix[ell[z3e]]].irac1flux), alog10(objectnew[ix[ell[z3e]]].irac4flux / objectnew[ix[ell[z3e]]].irac2flux), psym = 8, symsize=0.3,color = redcolor
;
oplot, alog10(objectnew[ix[ell[z4e]]].irac3flux / objectnew[ix[ell[z4e]]].irac1flux), alog10(objectnew[ix[ell[z4e]]].irac4flux / objectnew[ix[ell[z4e]]].irac2flux), psym = 8, symsize=0.3,color = bluecolor
;
oplot, alog10(objectnew[ix[ell[z5e]]].irac3flux / objectnew[ix[ell[z5e]]].irac1flux), alog10(objectnew[ix[ell[z5e]]].irac4flux / objectnew[ix[ell[z5e]]].irac2flux), psym = 8, symsize=0.3,color = purplecolor
;

oplot, alog10(objectnew[ix[star]].irac3flux / objectnew[ix[star]].irac1flux), alog10(objectnew[ix[star]].irac4flux / objectnew[ix[star]].irac2flux), psym = 8,  symsize=0.3

xyouts, -1.4, 0.9, ". 0.1 < z < 0.2", color = greencolor
xyouts, -1.4, 0.8, ". 0.2 < z > 0.4", color = yellowcolor
xyouts, -1.4, 0.7, ". 0.4 < z > 0.6", color = redcolor
xyouts, -1.4, 0.6, ". 0.6 < z > 0.8", color = bluecolor
xyouts, -1.4, 0.5, ". 0.8 < z > 1.0", color = purplecolor

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Ell2_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = redcolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = redcolor

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Ell13_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = bluecolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = bluecolor

;---------------------

plot, alog10(objectnew[ix[sp[z1s]]].irac3flux / objectnew[ix[sp[z1s]]].irac1flux), alog10(objectnew[ix[sp[z1s]]].irac4flux / objectnew[ix[sp[z1s]]].irac2flux), psym = 8, xrange = [-1.5, 1.], yrange=[-1.0, 1.0],/xst,/yst, xtitle='log(S5.8/S3.6)', ytitle = 'log(S8.0/S4.5)', symsize=0.3 ,/nodata, title = "spirals as a function of z"; ticklen = 1,

oplot, alog10(objectnew[ix[sp[z1s]]].irac3flux / objectnew[ix[sp[z1s]]].irac1flux), alog10(objectnew[ix[sp[z1s]]].irac4flux / objectnew[ix[sp[z1s]]].irac2flux), psym = 8, symsize=0.3, color=greencolor

oplot, alog10(objectnew[ix[sp[z2s]]].irac3flux / objectnew[ix[sp[z2s]]].irac1flux), alog10(objectnew[ix[sp[z2s]]].irac4flux / objectnew[ix[sp[z2s]]].irac2flux), psym = 8, symsize=0.3,color = yellowcolor

oplot, alog10(objectnew[ix[sp[z3s]]].irac3flux / objectnew[ix[sp[z3s]]].irac1flux), alog10(objectnew[ix[sp[z3s]]].irac4flux / objectnew[ix[sp[z3s]]].irac2flux), psym = 8, symsize=0.3,color = redcolor

oplot, alog10(objectnew[ix[sp[z4s]]].irac3flux / objectnew[ix[sp[z4s]]].irac1flux), alog10(objectnew[ix[sp[z4s]]].irac4flux / objectnew[ix[sp[z4s]]].irac2flux), psym = 8, symsize=0.3,color = bluecolor

oplot, alog10(objectnew[ix[sp[z5s]]].irac3flux / objectnew[ix[sp[z5s]]].irac1flux), alog10(objectnew[ix[sp[z5s]]].irac4flux / objectnew[ix[sp[z5s]]].irac2flux), psym = 8, symsize=0.3,color = purplecolor

oplot, alog10(objectnew[ix[star]].irac3flux / objectnew[ix[star]].irac1flux), alog10(objectnew[ix[star]].irac4flux / objectnew[ix[star]].irac2flux), psym = 8,  symsize=0.3


xyouts, -1.4, 0.9, ". 0.1 < z < 0.2", color = greencolor
xyouts, -1.4, 0.8, ". 0.2 < z > 0.4", color = yellowcolor
xyouts, -1.4, 0.7, ". 0.4 < z > 0.6", color = redcolor
xyouts, -1.4, 0.6, ". 0.6 < z > 0.8", color = bluecolor
xyouts, -1.4, 0.5, ". 0.8 < z > 1.0", color = purplecolor
xyouts, -1.4, 0.4, "stars"
;plothist, objectnew.zphot, xrange=[0,2] , bin=0.05

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Sa_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = redcolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = redcolor

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Sd_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = bluecolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = bluecolor
;---------------------

r = where(alog10(objectnew[ix].irac4flux  / objectnew[ix].irac2flux) gt 1.5)
print, 'r', n_elements(r)
;plothyperz, ix[r]



ps_close, /noprint,/noid

;plothyperz, ix[agn[agninterest]]

end



; plot, objectnew[ix].irac3mag - objectnew[ix].irac4mag, objectnew[ix].irac1mag - objectnew[ix].irac2mag, psym = 8, xrange = [-1.2, 2.9], yrange=[-0.75, 1.0],/xst,/yst,symsize=0.4, xtitle='ch3 - ch4', ytitle = 'ch1 - ch2'
; oplot,[1,1]*(-0.1),[-0.2,1]
; oplot,[-0.1,1.0],[-0.2,0]
; oplot,[1.0,1.3],[0.0,1.0]
; oplot,[1,1]*(-0.1),[-0.2,1]
;
;
;plot, objectnew[ix].irac3mag - objectnew[ix].irac1mag, objectnew[ix].irac4mag - objectnew[ix].irac2mag, psym = 8, xrange = [-1.2, 2.;9], yrange=[-1.0, 2.0],/xst,/yst,symsize=0.4, xtitle='ch3 - ch1', ytitle = 'ch4 - ch2'


;plothist, objectnew.rfwhm, bin = 0.1,xrange=[1,20], /noprint
;targ = where(objectnew.acsmag lt 23 and objectnew.acsmag gt 0 )
;print, "n_elements(targ)", n_elements(targ)

;plot, objectnew[targ].acsfwhm, objectnew[targ].acsellip, xrange=[1,20], psym = 2, xtitle='acs fwhm', ytitle = 'acs ellipticity'

;jx = where(objectnew.irac1mag lt 90.0 and objectnew.irac1mag gt 0.0 and objectnew.irac4mag gt 0.0 and objectnew.irac4mag lt 90. and objectnew.irac3mag gt 0.0 and objectnew.irac3mag lt 90. and objectnew.irac2mag gt 0.0 and objectnew.irac2mag lt 90. and objectnew.acsfwhm le 3.0 and objectnew.acsellip lt 0.15)
;print, "n_elenents(jx)", n_elements(jx)

;plot, objectnew[jx].irac3mag - objectnew[jx].irac4mag, objectnew[jx].irac1mag - objectnew[jx].irac2mag, psym = 8, xrange = [-1.2, 2.9], yrange=[-0.75, 1.0],/xst,/yst,symsize=0.4, xtitle='ch3 - ch4', ytitle = 'ch1 - ch2'
; oplot,[1,1]*(-0.1),[-0.2,1]
; oplot,[-0.1,1.0],[-0.2,0]
; oplot,[1.0,1.3],[0.0,1.0]
; oplot,[1,1]*(-0.1),[-0.2,1]
