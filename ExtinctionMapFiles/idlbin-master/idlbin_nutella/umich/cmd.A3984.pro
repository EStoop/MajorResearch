PRO cmd
close,/all
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; february 2004
;
;this code makes a color magnitude diagram of "galaxies" in A3984
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

numoobjects = 45000
rgalaxy = replicate({object, xcenter:0D,ycenter:0D,fluxa:0D,maga:0D,fluxb:0D,magb:0D,fwhm:0D,isoarea:0D},numoobjects)
Vgalaxy = replicate({Vobject, Vxcenter:0D,Vycenter:0D,Vfluxa:0D,Vmaga:0D,Vfluxb:0D,Vmagb:0D,Vfwhm:0D,Visoarea:0D},numoobjects)
mgalaxy = replicate({what, xcen:0D,ycen:0D, vr:0D,rmag:0D, member:0D},numoobjects)

i = 0
rsum = 0.
Vsum = 0.

openr, lun, "/n/Godiva4/jkrick/A3984/cmd/SExtractor.r.cat", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF,lun,junk, rx, ry, fluxaper,magaper,erraper,fluxbest,magbest,errbest,fwhm,rra,rdec,iso
    IF (fwhm GT 4.5) AND (fluxbest GT 0.0) AND (iso GT 6.0) THEN BEGIN   
      ;no bother including the stars and junk
        rgalaxy[i] = {object,  rx, ry, fluxaper,magaper,fluxbest,magbest,fwhm,iso}
        i = i + 1
    ENDIF
ENDWHILE

rgalaxy = rgalaxy[0:i - 1]
print, "there are i,",i," r galaxies"
close, lun
free_lun, lun


;read in the V magnitudes and match them up to make a V-r

m = 0
ship = 0
nship = 0
vra = 0D
vdec = 0D
xshift = 0.;-190
yshift = 0.;+60
openr, lun2, "/n/Godiva4/jkrick/A3984/cmd/SExtractor.V.cat", /GET_LUN

WHILE (NOT eof(lun2)) DO BEGIN
    READF, lun2,junk, vx, vy, vfluxaper,vmagaper,verraper,vfluxbest,vmagbest,$
     verrbest,vfwhm, vra,vdec,viso
;    xyad,header,vx,vy,vra,vdec      ;find the ra and dec from the better astrometry image
    membership = 0.                  ;assume it is not a member

    IF (vfwhm GT 4.2) AND (vfluxbest GT 0.0) AND (viso GT 6.0)THEN BEGIN ;not star or junk
        Vgalaxy[i] = {Vobject,  vx, vy, vfluxaper,vmagaper,vfluxbest,vmagbest,vfwhm,viso}
        FOR j = 0, i - 1, 1 DO BEGIN


            IF (vx LT rgalaxy[j].xcenter + 0.6 +xshift) AND $
              (vx GT rgalaxy[j].xcenter - 0.6 + xshift) AND $
              (vy LT rgalaxy[j].ycenter + 0.6 + yshift) AND $ 
              (vy GT rgalaxy[j].ycenter - 0.6 + yshift)  THEN BEGIN
                
                mgalaxy[m] = {what, vx,vy,vmagbest - rgalaxy[j].magb ,rgalaxy[j].magb,membership}

                IF rgalaxy[j].magb LT 19.0  AND vmagbest - rgalaxy[j].magb GT 0.4 THEN  $
                    print, "brightblue", vx,vy, rgalaxy[j].magb, vmagbest - rgalaxy[j].magb
;                IF rgalaxy[j].magb LT 23.0  AND vmagbest - rgalaxy[j].magb GT 1.0 THEN  $
;                    print, "brightblue", vx,vy, rgalaxy[j].magb, vmagbest - rgalaxy[j].magb

                m = m + 1
            ENDIF

        ENDFOR
    ENDIF
    
ENDWHILE
mgalaxy = mgalaxy[0:m-1]
close,lun2
free_lun, lun2

print, "There are ",m," objects in this plot"


;-------------------------------------------------------------------
;fitting section
;---------------------------------------------------------------------

;sort them in rmag space to only fit to the better data
sortindex= sort(mgalaxy.rmag)
rmagsort= mgalaxy.rmag[sortindex]
vrsort = mgalaxy.vr[sortindex]

;print, mgalaxy.rmag

r1 = fltarr(m)
vr1 = fltarr(m)
r2 = fltarr(m)
vr2 = fltarr(m)
t= 0
u = 0
;make the cutoff at r = 21, fit before 21 and after 21 with two
;different functions
;also do not fit anything brightwards of 16.1
FOR se = 0, m-1, 1 DO begin
    IF (mgalaxy[se].rmag LE 24) AND (mgalaxy[se].rmag GT 17.5) THEN BEGIN
        r1(t) = mgalaxy[se].rmag
        vr1(t) = mgalaxy[se].vr
        t = t +1
    ENDIF 
    IF (mgalaxy[se].rmag LT 21) AND (mgalaxy[se].rmag GT 17.5) THEN BEGIN
        r2(u) = mgalaxy[se].rmag
        vr2(u) = mgalaxy[se].vr
        u = u + 1
    ENDIF

ENDFOR
print, t, " the number of galaxies < 21"
r1 = r1[0:t - 1]
vr1 = vr1[0:t - 1]
r2 = r2[0:u-1]
vr2 = vr2[0:u-1]

;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( mgalaxy.rmag, mgalaxy.vr, yfit, sig, coeff_sig)
coeff1 = ROBUST_LINEFIT(r1,vr1, yfit1, sig1, coeff_sig1)
coeff2 = ROBUST_LINEFIT(r2,vr2, yfit2, sig2, coeff_sig2)

err = dindgen(40) - dindgen(40) + 1
start = [-0.02,1.0]
r1j = fltarr(40)
y1j = fltarr(40)
r1j = r1[3:42]
y1j = yfit1[3:42]
sortindex = sort(r1j)
sortr1 = r1j[sortindex]
sortyfit1 = y1j[sortindex]
result = MPFITFUN('linear',sortr1, sortyfit1,err, start)
print, "result", result

;--------------------------------------------------------
;for every galaxy , given it's r value, does it's v-r value lie
;within the acceptable range?
;-----------------------------------------------------------
openw, nolun, "/n/Godiva4/jkrick/A3984/cmd/member.txt",/get_lun
openw, newlun, "/n/Godiva4/jkrick/A3984/cmd/density.txt",/get_lun
memgalincut = 0.
nmemgalincut = 0.
memgalnincut = 0.
nmemgalnincut = 0.
galcount = 0.
num = 0
memarray= fltarr(20000)
memcount = 0
density = 0
xarr = fltarr(20000)
yarr = fltarr(20000)

FOR a = 0, m-1 ,1 DO BEGIN

    IF (mgalaxy[a].member GT 0 ) AND ( mgalaxy[a].member LT 9.)THEN num = num + 1

     ;is v-r LT yfit[a] + sig and GT yfit[a] - sig
    limithigh = ((result(0))*mgalaxy[a].rmag) + result(1)+ sig2 
    limitlow = ((result(0))*mgalaxy[a].rmag) + result(1)- sig2 
    IF (mgalaxy[a].vr LT limithigh )AND (mgalaxy[a].vr GT limitlow) THEN BEGIN
       ;galaxy "a" makes the color cut
        galcount = galcount + 1
        printf, nolun, mgalaxy[a].xcen,mgalaxy[a].ycen, mgalaxy[a].rmag, $
          mgalaxy[a].vr + mgalaxy[a].rmag

        memarray[memcount] = mgalaxy[a].rmag
        xarr[memcount] = mgalaxy[a].xcen
        yarr[memcount] = mgalaxy[a].ycen

        memcount = memcount + 1
        ;where are those central galaxies?
;        IF mgalaxy[a].xcen LT 1600 AND mgalaxy[a].xcen GT 1250 AND $
;        mgalaxy[a].ycen LT 1680 AND mgalaxy[a].ycen GT 1230 THEN
;        print, $
;        "galaxy makes cut", mgalaxy[a].xcen, mgalaxy[a].ycen
        
        IF (mgalaxy[a].member GT 0.) AND ( mgalaxy[a].member LT 9.)THEN BEGIN
            memgalincut = memgalincut + 1
        ENDIF
        IF (mgalaxy[a].member GT 9.) THEN BEGIN
            nmemgalincut = nmemgalincut + 1
        ENDIF
        dist1 = sqrt((mgalaxy[a].xcen-1161.9)^2 + (mgalaxy[a].ycen - 2150.)^2) 
;        dist2 = sqrt((mgalaxy[a].xcen-1157.)^2 + (mgalaxy[a].ycen - 2267.)^2) 

        IF dist1 LE 1097. THEN BEGIN
            rsum = rsum + 10^((24.6 - mgalaxy[a].rmag)/2.5)
            Vsum = Vsum + 10^((24.3 - (mgalaxy[a].vr + mgalaxy[a].rmag))/2.5)
        ENDIF
        IF dist1 LE 1011 AND mgalaxy[a].rmag LT 21.17 THEN BEGIN
            density = density + 1
            printf, newlun, mgalaxy[a].xcen, mgalaxy[a].ycen, mgalaxy[a].rmag
        ENDIF
        
;        ENDIF

    ENDIF ELSE BEGIN
        IF (mgalaxy[a].member GT 0.) AND ( mgalaxy[a].member LT 9.)THEN BEGIN
            memgalnincut = memgalnincut + 1
        ENDIF
         IF (mgalaxy[a].member GT 9.) THEN BEGIN
            nmemgalnincut = nmemgalnincut + 1
        ENDIF
        ;galaxy"a" does not make the color cut
        ;will want to mask these galaxies to
        ;determine a cluster light profile
    ENDELSE
ENDFOR
memarray = memarray[0:memcount -1]
xarr = xarr[0:memcount -1]
yarr = yarr[0:memcount -1]

sortarray = sort(memarray)
memarray = memarray[sortarray]
xarr = xarr[sortarray]
yarr = yarr[sortarray]
print, "sorted magnitudes", memarray[0],xarr[0],yarr[0],memarray[1],xarr[1],yarr[1],memarray[2],xarr[2],yarr[2],memarray[3],xarr[3],yarr[3],$
memarray[4],xarr[4],yarr[4],memarray[5],xarr[5],yarr[5]

;density = 0
;FOR num = 0, n_elements(memarray) - 1, 1 DO BEGIN
;    dist1 = sqrt((xarr(num)-1162.)^2 + (yarr(num) - 2150.)^2) 
;    IF (memarray(num) LT ( memarray[9] + 3)) AND (dist1 LE 1056.) THEN density = density + 1
;ENDFOR
;-------------------------------------------------------------
;-------------------------------------------------------------
;for plotting up the 2 different clumps
;----------------------------------------------------------------
newrarr = fltarr(20000)
newvrarr = fltarr(20000)
newnrarr = fltarr(20000)
newnvrarr = fltarr(20000)
p = 0
q = 0
FOR g = 0, m - 1, 1 DO BEGIN
;    IF (sqrt((mgalaxy[g].xcen - 1605.)^2 + (mgalaxy[g].ycen - 1218)^2) LT 198 ) AND (mgalaxy[g].rmag LT 21 ) THEN BEGIN
    IF (mgalaxy[g].ycen GT 1730) AND (mgalaxy[g].rmag LT 21 ) THEN BEGIN
        ;nearer to the cD clump
        newrarr[p] = mgalaxy[g].rmag
        newvrarr[p] = mgalaxy[g].vr
        p = p+1
    ENDIF
;    IF (sqrt((mgalaxy[g].xcen - 1322.)^2 + (mgalaxy[g].ycen - 1470)^2) LT 182 ) AND (mgalaxy[g].rmag LT 21 ) THEN BEGIN
    IF (mgalaxy[g].ycen LT 1730) AND (mgalaxy[g].rmag LT 21 ) THEN BEGIN
        ;nearer to the bottom clump
        newnrarr[q] = mgalaxy[g].rmag
        newnvrarr[q] = mgalaxy[g].vr
        q = q+1
    ENDIF
    
        
ENDFOR

newrarr =   newrarr[0:p-1]
newvrarr = newvrarr[0:p-1]
newnrarr = newnrarr[0:q-1]
newnvrarr = newnvrarr[0:q-1]
;#######################################


device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

;device, filename = '/n/Godiva1/jkrick/A3888/cmd/cmdnocolor.ps', /portrait, $
;  BITS=8, scale_factor=0.9 , /color
device, filename = '/n/Godiva4/jkrick/A3984/cmd/junk.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot,[0,0,0],[0,0,0],  thick = 3, psym = 2, xthick = 3, ythick = 3,charthick = 3,$
                        xrange = [15,27], yrange = [-1,1.5], ystyle = 1,xstyle = 1,$
	  ytitle = "V-r", $;title = 'cmd looking for cluster galaxies', $
                        xtitle = "r (mag)", charsize = 1.5
oplot,mgalaxy.rmag, mgalaxy.vr,  thick = 3, psym = 2,color = colors.darkgray
;;oplot, newrarr, newvrarr, psym = 5, color = colors.black, thick = 3
;;oplot, newnrarr,newnvrarr, psym = 6 ,color = colors.black, thick = 3

;oplot,r1, yfit1, thick = 3;, color = colors.yellow
;oplot, r1, yfit1 + sig1;, thick = 3, color = colors.yellow
;oplot, r1, yfit1 - sig1;, thick = 3, color = colors.yellow

xvals = findgen(11) +16
oplot, xvals, ((result(0))*xvals) + result(1), thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)+ sig2 , thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)- sig2 , thick = 3;, color = colors.orange

;plothist, mgalaxy.vr, bin = 0.05, peak = 1.0,  thick = 3, xthick = 3, ythick = 3,charthick = 3,$
;	   title = 'V-r colors of galaxies in A3984 field', $
;                        xtitle = "V -r", xrange = [-1,1]
;plothist, newnvrarr, bin = 0.05, peak = 1.0,  thick = 3, xthick = 3, ythick = 3,charthick = 3,$
;	   title = 'V-r colors of galaxies in A3984 field', $
;                        xtitle = "V -r", color = colors.blue, xrange = [-1,1]
;plothist, newvrarr,bin = 0.05, peak = 1.0,  thick = 3, xthick = 3, ythick = 3,charthick = 3,$
;	   title = 'V-r colors of galaxies in A3984 field', $
;                        xtitle = "V -r", color = colors.red, xrange = [-1,1]

;phist, newnvrarr, -0.5,1.0,0.01,/cumul

;phist, newvrarr, -0.5,1.0,0.01, /cumul,/overplot

print, newvrarr

kstwo, newnvrarr,newvrarr, D, prob   
print,"D,prob",D,prob
kstwo, newnvrarr,mgalaxy.vr, D, prob   
print,"D,prob",D,prob

print, "sig2", sig2
;plothist,mgalaxy.vr,xhist,yhist,bin=0.05
print,"mean", mean(mgalaxy.vr)
device, /close
set_plot, mydevice

undefine, rgalaxy
;undefine, lgalaxy
undefine, mgalaxy
undefine, Vgalaxy


print, rsum, Vsum
print, "density", density
close, nolun
free_lun, nolun
close, newlun
free_lun, newlun
END


;           galaxy[m].vr= vmagbest - rgalaxy[j].magb 
;            galaxy[m].rmag = rgalaxy[j].magb 
;            galaxy[m].xcen = vx
;            galaxy[m].ycen = vy
;xyouts, 16.5,-0.6,"confirmed members", charthick = 5, charsize = 1.2
;xyouts, 16.5,-0.7, "confirmed non-members", charthick = 5, charsize = 1.2


;plothist, mgalaxy.vr, bin = 0.05, peak = 1.0,  thick = 3, xthick = 3, ythick = 3,charthick = 3,$
;	   title = 'V-r colors of galaxies in A3888 field', $
;                        xtitle = "V -r"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;determine how many galaxies there are per magnitude bin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;sort Vgalaxy on magnitude
;sortindex = sort(Vgalaxy.vmagb)
;sortedmag = Vgalaxy.vmagb[sortindex]

;plothist, sortedmag, xhist, yhist, bin = 0.5
;decide on bin size = 0.5mags
;numogals = 0
;bin = fltarr(19)
;mag = fltarr(19)
;j = 0
;for b = 15, 24, 0.5 DO BEGIN               ;count number of galaxies
;WHILE (sortedmag(j) LT b + 0.5)  DO BEGIN
;    numogals = numogals + 1
;    print, "numogals + 1", sortedmag(j)
;    j = j + 1
;ENDWHILE
;bin(i) = 

