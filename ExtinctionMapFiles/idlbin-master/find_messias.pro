pro find_messias
; put the messias wedges onto the plots
; get selection fractions of each of these selections and also lacy and stern.
;but lacy and stern use 3.6 and 5.8.
;is wirck in AB?

;ps_start, filename= '/Volumes/irac_drive/nutella/messias/agn_selection.ps'
ps_start, filename= '~/nutella/messias/agn_selection.ps'
!P.multi = [0,0,1]
!P.charthick=3
!P.thick=3
!X.thick=3
!Y.thick=3

vsym, /polygon, /fill

;close, /all
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("green", !D.Table_Size-4)

;!P.multi = [0,3,1]
restore, '/Users/jkrick/nutella/idlbin/objectnew_akari.sav'
;restore, '/Volumes/irac_drive/nutella/idlbin/objectnew.sav'

KI = where(objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 );and (objectnew.wirckmag - objectnew.irac2mag) gt 0 and (objectnew.irac2mag - objectnew.irac4mag) gt 0)

print, n_elements(KI)

;need to convert kmag to AB magnitudes
wirckab = 666.8*10^(objectnew[KI].wirckmag/(-2.5))
objectnew[KI].wirckmag = -2.5*alog10(wirckab) +8.926


plot, objectnew[KI].irac2mag - objectnew[KI].irac4mag, objectnew[KI].wirckmag - objectnew[KI].irac2mag,  psym = 8, xtitle = '4.5 - 8.0', ytitle  = 'K-4.5', title = string(n_elements(KI)) + ' objects detected at K, ch2, & ch4', charsize = 1, xrange = [-2,2], yrange = [-2,2];, symsize = 0.3

;now ask which ones does hyperz think are AGN
KIagn = where(objectnew[KI].spt eq 6 or objectnew[KI].spt eq 7 or objectnew[KI].spt eq 13$
or objectnew[KI].spt eq 14 or objectnew[KI].spt eq 15)

print, 'n kiang', n_elements(KIagn)

KInoagn = where(objectnew[KI].spt eq 1 or objectnew[KI].spt eq 2 or objectnew[KI].spt eq 3 or objectnew[KI].spt eq 4 or objectnew[KI].spt eq 5 or objectnew[KI].spt eq 8 or objectnew[KI].spt eq 9 or objectnew[KI].spt eq 10 or objectnew[KI].spt eq 11 or objectnew[KI].spt eq 12)

oplot, objectnew[KI[kiagn]].irac2mag - objectnew[KI[kiagn]].irac4mag, objectnew[KI[kiagn]].wirckmag - objectnew[KI[kiagn]].irac2mag, psym = 8,color = bluecolor;, symsize = 0.3

oplot, [0,0], [0,2.0], linestyle = 2, thick =3
oplot, [0,2.0], [0, 0], linestyle = 2, thick =3



;;-----------------------

KIM = where(objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.mips24flux gt 17.3 and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 );and objectnew.zphot gt 3.0 and objectnew.zphot lt 6.0 );and (objectnew.wirckmag - objectnew.irac2mag) gt 0 and (objectnew.irac2mag - objectnew.irac4mag) gt 0)

print, n_elements(KIM)

plot, objectnew[KIM].irac2mag - objectnew[KIM].irac4mag, objectnew[KIM].wirckmag - objectnew[KIM].irac2mag, psym = 8, xtitle = '4.5 - 8.0', ytitle  = 'K-4.5', title = string(n_elements(KIM)) + ' objects detected at K, ch2, ch4, & 24', charsize = 1,xrange = [-2,2], yrange = [-2,2];, symsize = 0.3
;now as which ones does hyperz think are AGN
KIMagn = where(objectnew[KIM].spt eq 6 or objectnew[KIM].spt eq 7 or objectnew[KIM].spt eq 13$
or objectnew[KIM].spt eq 14 or objectnew[KIM].spt eq 15)

KIMnoagn = where(objectnew[KIM].spt eq 1 or objectnew[KIM].spt eq 2 or objectnew[KIM].spt eq 3 or objectnew[KIM].spt eq 4 or objectnew[KIM].spt eq 5 or objectnew[KIM].spt eq 8 or objectnew[KIM].spt eq 9 or objectnew[KIM].spt eq 10 or objectnew[KIM].spt eq 11 or objectnew[KIM].spt eq 12)

oplot,  objectnew[KIM[kimagn]].irac2mag - objectnew[KIM[kimagn]].irac4mag, objectnew[KIM[kimagn]].wirckmag - objectnew[KIM[kimagn]].irac2mag, psym = 8,color = bluecolor;, symsize = 0.3

oplot, [0,0], [0,2.0], linestyle = 2, thick =3
oplot, [0,2.0], [0, 0], linestyle = 2, thick =3

x = where(objectnew[kim].xnetcounts gt 0)
oplot, objectnew[KIM[x]].irac2mag - objectnew[KIM[x]].irac4mag, objectnew[KIM[x]].wirckmag - objectnew[KIM[x]].irac2mag, psym = 8,color = redcolor, symsize = 0.5;


;;-----------------------


mips24abmag_KIM = 23.2 - 2.5*alog10(objectnew[KIM].mips24flux)

plot, objectnew[KIM].irac2mag - objectnew[KIM].irac4mag, objectnew[KIM].irac4mag - mips24abmag_KIM, psym  = 8, xtitle = '4.5 - 8.0', ytitle = '8.0 - 24', charsize = 1, title = string(n_elements(KIM)) + ' objects detected at K, ch2, ch4, & 24', yrange = [-2, 4], xrange = [-1.5,3.0];, symsize = 0.3

oplot, objectnew[KIM[kimagn]].irac2mag - objectnew[KIM[kimagn]].irac4mag, objectnew[KIM[kimagn]].irac4mag - mips24abmag_KIM[kimagn], psym  = 8, color = bluecolor;, symsize = 0.3

;oplot, [-2.0,2.0], [0.5, 0.5], linestyle = 2, thick =3
;x = [-2,2]
;y =[-4,4]
;oplot, (y - 2.465) / (-3.275), linestyle = 2, thick = 3

oplot, [0.6,3.0], [0.5, 0.5], linestyle = 2, thick = 3
oplot, [-0.6, 0.6], [4.1, 0.5], linestyle = 2, thick = 3



oplot, objectnew[KIM[x]].irac2mag - objectnew[KIM[x]].irac4mag, objectnew[KIM[x]].irac4mag - mips24abmag_KIM[x], psym  = 8, color = redcolor, symsize = 0.5

n = where(objectnew[KIM].wirckmag - objectnew[KIM].irac2mag gt 0)
oplot, objectnew[KIM[n]].irac2mag - objectnew[KIM[n]].irac4mag, objectnew[KIM[n]].irac4mag - mips24abmag_KIM[n], psym  = 6, color = greencolor

;-----------------KI
;now how many objects are in each of the selections?
; what fraction of those are AGN?
; what fraction of AGN are outside of the selection?

KI_messias = where(objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and (objectnew.wirckmag - objectnew.irac2mag) gt 0 and (objectnew.irac2mag - objectnew.irac4mag) gt 0)


KI_messias_noagn = where(objectnew[KI_messias].spt eq 1 or objectnew[KI_messias].spt eq 2 or objectnew[KI_messias].spt eq 3 or objectnew[KI_messias].spt eq 4 or objectnew[KI_messias].spt eq 5 or objectnew[KI_messias].spt eq 8 or objectnew[KI_messias].spt eq 9 or objectnew[KI_messias].spt eq 10 or objectnew[KI_messias].spt eq 11 or objectnew[KI_messias].spt eq 12)

; how many total AGN are there

print, 'fraction of AGN which are selected by KI', n_elements(KI_messias), n_elements(KIagn)
print, 'fraction of those selected which are not AGN', n_elements(KI_messias_noagn), n_elements(KI_messias)


;-----------------KIM
;now how many objects are in each of the selections?
; what fraction of those are AGN?
; what fraction of AGN are outside of the selection?

KIM_messias = where(objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.mips24flux gt 17.3 and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and (objectnew.wirckmag - objectnew.irac2mag) gt 0 and (objectnew.irac4mag - (23.2 - 2.5*alog10(objectnew.mips24flux))) gt 0.5 and (objectnew.irac4mag - (23.2 - 2.5*alog10(objectnew.mips24flux))) gt (-3.275*(objectnew.irac2mag - objectnew.irac4mag) + 2.465))

print, 'KIM_messias', n_elements(KIM_messias)
KIM_messias_noagn = where(objectnew[KIM_messias].spt eq 1 or objectnew[KIM_messias].spt eq 2 or objectnew[KIM_messias].spt eq 3 or objectnew[KIM_messias].spt eq 4 or objectnew[KIM_messias].spt eq 5 or objectnew[KIM_messias].spt eq 8 or objectnew[KIM_messias].spt eq 9 or objectnew[KIM_messias].spt eq 10 or objectnew[KIM_messias].spt eq 11 or objectnew[KIM_messias].spt eq 12)

; how many total AGN are there

print, 'fraction of AGN which are selected by KIM', n_elements(KIM_messias), n_elements(KIMagn)
print, 'fraction of those selected which are not AGN', n_elements(KIM_messias_noagn), n_elements(KIM_messias)



;-------------------------------------------------
;now look at the akari photometry
a = where(objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and objectnew.flux11 gt 0 and objectnew.flux18 gt 0)
print, 'n_akari', n_elements(a)

;convert flux to ab mag
magab11 = 23.2 - 2.5*alog10(objectnew[a].flux11)
magab18 = 23.2 - 2.5*alog10(objectnew[a].flux18)

plot, objectnew[a].irac2mag - magab11, magab11 - magab18, psym = 8, xtitle = '4.5 - 11', ytitle = '11 - 18',  title = string(n_elements(a)) + ' objects detected at 4.5, 11, & 18'

oplot, [0.4, 5.0], [0.3,0.3] ,linestyle = 2, thick = 3
oplot, [-0.5, 0.4], [4.0,0.3] ,linestyle = 2, thick = 3

aagn = where(objectnew[a].spt eq 6 or objectnew[a].spt eq 7 or objectnew[a].spt eq 13$
or objectnew[a].spt eq 14 or objectnew[a].spt eq 15)

print, 'n aagn', n_elements(aagn)

anoagn = where(objectnew[a].spt eq 1 or objectnew[a].spt eq 2 or objectnew[a].spt eq 3 or objectnew[a].spt eq 4 or objectnew[a].spt eq 5 or objectnew[a].spt eq 8 or objectnew[a].spt eq 9 or objectnew[a].spt eq 10 or objectnew[a].spt eq 11 or objectnew[a].spt eq 12)

oplot, objectnew[a[aagn]].irac2mag - magab11[aagn], magab11[aagn] - magab18[aagn],psym = 8,color = bluecolor

x = where(objectnew[a].xnetcounts gt 0)
print, 'x', n_elements(x)
print, objectnew[a[x]].irac2mag - magab11[x]
print, 'y', magab11[x] - magab18[x]
oplot,objectnew[a[x]].irac2mag - magab11[x], magab11[x] - magab18[x], psym = 8,color = redcolor, symsize = 0.5;



ps_end

end
