PRO cmd
close,/all
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; february 2004
;
;this code makes a color magnitude diagram of "galaxies" in A3888
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fits_read,"/n/Godiva1/jkrick/A3888/cmd/V.wcs.fits",data,header
;fits_read,"/n/Godiva1/jkrick/A3888/fullV.wcs.fits",data,header

numoobjects = 45000
rgalaxy = replicate({object, xcenter:0D,ycenter:0D,fluxa:0D,maga:0D,fluxb:0D,magb:0D,fwhm:0D,isoarea:0D},numoobjects)
Vgalaxy = replicate({Vobject, Vxcenter:0D,Vycenter:0D,Vfluxa:0D,Vmaga:0D,Vfluxb:0D,Vmagb:0D,Vfwhm:0D,Visoarea:0D},numoobjects)
mgalaxy = replicate({what, xcen:0D,ycen:0D, vr:0D,rmag:0D, member:0D},numoobjects)

i = 0

openr, lun, "/n/Godiva1/jkrick/A3888/cmd/SExtractor2.r.cat", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF,lun,junk, rx, ry, fluxaper,magaper,erraper,fluxbest,magbest,errbest,fwhm,rra,rdec,iso
    IF (fwhm GT 5.1) AND (fluxbest GT 0.0) AND (iso GT 6.0) THEN BEGIN   ;no bother including the stars and junk
        rgalaxy[i] = {object,  rx, ry, fluxaper,magaper,fluxbest,magbest,fwhm,iso}
        i = i + 1
    ENDIF
ENDWHILE

rgalaxy = rgalaxy[0:i - 1]
print, "there are i,",i," r galaxies"
close, lun
free_lun, lun

;read in the LARCS membership information
numogals = 218
l = 0
lgalaxy = replicate({larcs,ra:0D,dec:0D},numogals)
openr, lun3, "/n/Godiva1/jkrick/A3888/pimbblet/larcs.cz.txt", /GET_LUN
WHILE (NOT eof(lun3)) DO BEGIN
    READF, lun3,junk, junk,lra,ldec,junk,junk
    lgalaxy[l] = {larcs,lra,ldec}
    l = l +1
ENDWHILE
close, lun3
free_lun, lun3

;read in the LARCS non-membership information
numogals = 208
n = 0
ngalaxy = replicate({nlarcs,ra:0D,dec:0D},numogals)
openr, lun4, "/n/Godiva1/jkrick/A3888/pimbblet/nonmembers.txt", /GET_LUN
WHILE (NOT eof(lun4)) DO BEGIN
    READF, lun4,junk, junk,nra,ndec,junk,junk
    ngalaxy[n] = {larcs,nra,ndec}
    n = n +1
ENDWHILE
close, lun4
free_lun, lun4

;read in the teague membership information
numogals = 73
t = 0
tgalaxy = replicate({tcg,ra:0D,dec:0D},numogals)
openr, lun6, "/n/Godiva1/jkrick/A3888/cmd/teaguemembers.txt", /GET_LUN
WHILE (NOT eof(lun6)) DO BEGIN
    READF, lun6,tra,tdec
    tgalaxy[t] = {tcg,tra,tdec}
    t = t+1
ENDWHILE
close, lun6
free_lun, lun6

;read in the teague nonmembership information
numogals = 23
tn = 0
tngalaxy = replicate({tcgn,ra:0D,dec:0D},numogals)
openr, lun7, "/n/Godiva1/jkrick/A3888/cmd/teaguenonmem.txt", /GET_LUN
WHILE (NOT eof(lun7)) DO BEGIN
    READF, lun7,tnra,tndec
    tngalaxy[tn] = {tcgn,tnra,tndec}
    tn = tn+1
ENDWHILE
close, lun7
free_lun, lun7


;read in the V magnitudes and match them up to make a V-r

m = 0
ship = 0
nship = 0
vra = 0D
vdec = 0D
openr, lun2, "/n/Godiva1/jkrick/A3888/cmd/SExtractor2.V.cat", /GET_LUN
WHILE (NOT eof(lun2)) DO BEGIN
    READF, lun2,junk, vx, vy, vfluxaper,vmagaper,verraper,vfluxbest,vmagbest,verrbest,vfwhm, vra,vdec,viso
    xyad,header,vx,vy,vra,vdec                           ;find the ra and dec from the better astrometry image
    membership = 0.                                                ;assume it is not a member

    IF (vfwhm GT 5.1) AND (vfluxbest GT 0.0) AND (viso GT 6.0)THEN BEGIN ;not star or junk
        Vgalaxy[i] = {Vobject,  vx, vy, vfluxaper,vmagaper,vfluxbest,vmagbest,vfwhm,viso}
        FOR j = 0, i - 1, 1 DO BEGIN
            
            IF j LT l  THEN BEGIN                             ;is this v galaxy really a member according to LARCS?            
                IF (vra LT lgalaxy[j].ra + 0.001) AND (vra GT lgalaxy[j].ra - 0.001) AND $
                     (vdec LT lgalaxy[j].dec + 0.001) AND (vdec GT lgalaxy[j].dec - 0.001) THEN BEGIN
                    membership =  1.
                    ship = ship + 1
                ENDIF
            ENDIF

            IF j LT n  THEN BEGIN                             ;is this v galaxy really a nonmember according to LARCS?            
                IF (vra LT ngalaxy[j].ra + 0.001) AND (vra GT ngalaxy[j].ra - 0.001) AND $
                     (vdec LT ngalaxy[j].dec + 0.001) AND (vdec GT ngalaxy[j].dec - 0.001) THEN BEGIN
                    membership =  10.
                    nship = nship + 1
;                  print,vx,vy
                ENDIF
            ENDIF

            IF j LT t  THEN BEGIN ;is this v galaxy really a member according to TCG?            
                IF (vra LT tgalaxy[j].ra + 0.001) AND (vra GT tgalaxy[j].ra - 0.001) AND $
                     (vdec LT tgalaxy[j].dec + 0.001) AND (vdec GT tgalaxy[j].dec - 0.001) THEN BEGIN
                     ;but I don't want to recount it as a member if I already got it under the LARCS galaxies
                    IF (membership NE 1) THEN ship = ship + 1
                    membership =  1.
                ENDIF
            ENDIF

            IF j LT tn  THEN BEGIN                             ;is this v galaxy really a nonmember according to TCG?            
                IF (vra LT tngalaxy[j].ra + 0.001) AND (vra GT tngalaxy[j].ra - 0.001) AND $
                     (vdec LT tngalaxy[j].dec + 0.001) AND (vdec GT tngalaxy[j].dec - 0.001) THEN BEGIN
                    IF (membership NE 10) THEN nship = nship + 1      ;already got it with LARCS
                    membership =  10.
                ENDIF
            ENDIF


            IF (vx LT rgalaxy[j].xcenter + 0.6) AND  (vx GT rgalaxy[j].xcenter - 0.6) AND $
                 (vy LT rgalaxy[j].ycenter + 0.6) AND (vy GT rgalaxy[j].ycenter - 0.6)  THEN BEGIN
                mgalaxy[m] = {what, vx,vy,vmagbest - rgalaxy[j].magb ,rgalaxy[j].magb,membership}

;                IF rgalaxy[j].magb LT 20.0 AND vmagbest - rgalaxy[j].magb GT 0.6 THEN  $
;                     print, "red", vx,vy, rgalaxy[j].magb, vmagbest -
;                     rgalaxy[j].magb

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
;sortindex= sort(mgalaxy.rmag)
;rmagsort= mgalaxy.rmag[sortindex]
;vrsort = mgalaxy.vr[sortindex]

;print, mgalaxy.rmag

r1 = fltarr(m)
vr1 = fltarr(m)
r2 = fltarr(m)
vr2 = fltarr(m)
t= 0
u = 0
;make the cutoff at r = 21, fit before 21 and after 21 with two
;different functions
;also do not fit anything brightwards of 17.1
FOR se = 0, m-1, 1 DO begin
    IF (mgalaxy[se].rmag LE 21) AND (mgalaxy[se].rmag GT 17.1) THEN BEGIN
        r1(t) = mgalaxy[se].rmag
        vr1(t) = mgalaxy[se].vr
        t = t +1
    ENDIF ELSE BEGIN
        r2(u) = mgalaxy[se].rmag
        vr2(u) = mgalaxy[se].vr
        u = u + 1
    ENDELSE

ENDFOR
r1 = r1[0:t - 1]
vr1 = vr1[0:t - 1]
r2 = r2[0:u-1]
vr2 = vr2[0:u-1]

;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( mgalaxy.rmag, mgalaxy.vr, yfit, sig, coeff_sig)
coeff1 = ROBUST_LINEFIT(r1,vr1, yfit1, sig1, coeff_sig1)
coeff2 = ROBUST_LINEFIT(r2,vr2, yfit2, sig2, coeff_sig2)


;---------------------------------------------------------------
;did the galaxy make the color cut?
;-------------------------------------------------------------

;need to extrapolate yfit1 out to r = 25
;don't know what coeff means...
;so fit a line through some of the points
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


openw, nolun, "/n/Godiva1/jkrick/A3888/cmd/color2.txt",/get_lun
memgalincut = 0.
nmemgalincut = 0.
memgalnincut = 0.
nmemgalnincut = 0.
galcount = 0.
num = 0
density=0
densityg=0
rsum = 0.
Vsum = 0.
memarray= fltarr(20000)
memcount = 0
xarr = fltarr(20000)
yarr = fltarr(20000)
;for every galaxy , given it's r value, does it's v-r value lie
;within the acceptable range?

FOR a = 0, m-1 ,1 DO BEGIN

    IF (mgalaxy[a].member GT 0 ) AND ( mgalaxy[a].member LT 9.)THEN num = num + 1

     ;is v-r LT yfit[a] + sig and GT yfit[a] - sig
    limithigh = ((result(0))*mgalaxy[a].rmag) + result(1)+ sig1 
    limitlow = ((result(0))*mgalaxy[a].rmag) + result(1)- sig1 
    IF (mgalaxy[a].vr LT limithigh )AND (mgalaxy[a].vr GT limitlow) THEN BEGIN
       ;galaxy "a" makes the color cut
        galcount = galcount + 1
        memarray[memcount] = mgalaxy[a].rmag
        xarr[memcount] = mgalaxy[a].xcen
        yarr[memcount] = mgalaxy[a].ycen

        memcount = memcount + 1

        ;where are those central galaxies?
;        IF mgalaxy[a].xcen LT 1600 AND mgalaxy[a].xcen GT 1250 AND mgalaxy[a].ycen LT 1680 AND mgalaxy[a].ycen GT 1230 THEN print, "galaxy makes cut", mgalaxy[a].xcen, mgalaxy[a].ycen
        
        IF (mgalaxy[a].member GT 0.) AND ( mgalaxy[a].member LT 9.)THEN BEGIN
            memgalincut = memgalincut + 1
         ENDIF
         IF (mgalaxy[a].member GT 9.) THEN BEGIN
            nmemgalincut = nmemgalincut + 1
        ENDIF
       printf, nolun, mgalaxy[a].xcen,mgalaxy[a].ycen, mgalaxy[a].rmag
        dist1 = sqrt((mgalaxy[a].xcen- 1436.)^2 + (mgalaxy[a].ycen - 1524.)^2) 
        dist2 = sqrt((mgalaxy[a].xcen- 403.)^2 + (mgalaxy[a].ycen - 425.)^2) 

        IF dist1 LE 1345. THEN BEGIN
            rsum = rsum + 10^((24.6 - mgalaxy[a].rmag)/2.5)
            Vsum = Vsum + 10^((24.3 - (mgalaxy[a].vr + mgalaxy[a].rmag))/2.5)
        ENDIF

        IF dist1 LE 1182 AND mgalaxy[a].rmag LT 20.72 THEN density = density + 1;20.72
        IF dist2 LE 400 AND mgalaxy[a].rmag LT 20.72 THEN densityg = densityg + 1
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
;        printf, nolun, mgalaxy[a].xcen,mgalaxy[a].ycen, mgalaxy[a].rmag
    ENDELSE
ENDFOR

print, "percent of members included in color cut", memgalincut/(memgalincut+memgalnincut) *100.
print, "percent of nonmembers included in color cut", nmemgalincut/(nmemgalnincut + nmemgalincut) *100.
print, "percent of all objects that made the color cut", galcount / m *100.

print,"memgalincut", memgalincut
print,"memgalnincut", memgalnincut
print,"nmemgalincut", nmemgalincut
print,"nmemgalnincut", nmemgalnincut
print, "m = total number of galaxies", m
print,"galcount=total galaxies that made the cut", galcount
print, "total number of member galaxies", num

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
;    dist1 = sqrt((xarr(num)-1436.)^2 + (yarr(num) - 1524.)^2) 
;    IF (memarray(num) LT ( memarray[9] + 3)) AND (dist1 LE 1366.) THEN density = density + 1
;ENDFOR
;-------------------------------------------------------------
;for plotting up the confirmed members and non-members
;----------------------------------------------------------------
newrarr = fltarr(200)
newvrarr = fltarr(200)
newnrarr = fltarr(200)
newnvrarr = fltarr(200)
radarr = fltarr(m-1)
colorarr = fltarr(m-1)
count = 0
p = 0
q = 0
FOR g = 0, m - 1, 1 DO BEGIN
    IF (mgalaxy[g].member GT 0.) AND ( mgalaxy[g].member LT 9.)THEN BEGIN
        newrarr[p] = mgalaxy[g].rmag
        newvrarr[p] = mgalaxy[g].vr
        p = p+1
    ENDIF
    IF (mgalaxy[g].member GT 9.) THEN BEGIN
        newnrarr[q] = mgalaxy[g].rmag
        newnvrarr[q] = mgalaxy[g].vr
        q = q+1
    ENDIF
    
    ;get rid of saturated stars masquereding as bright blue galaxies
    IF (mgalaxy[g].rmag LT 16.0) THEN BEGIN
        mgalaxy[g].rmag = 0
        mgalaxy[g].vr = 0
    ENDIF
    IF (mgalaxy[g].rmag LT 16.5) AND (mgalaxy[g].vr  LT 0.0) THEN BEGIN
        mgalaxy[g].rmag = 0
        mgalaxy[g].vr = 0
    ENDIF

    ;make an array with all galaxies whose mag is LT 21
    IF (mgalaxy[g].rmag LT 22.) THEN BEGIN
        radarr[count] = sqrt((mgalaxy[g].xcen - 1556.)^2. + (mgalaxy[g].ycen - 1784.)^2.)
        colorarr[count] = mgalaxy[g].vr
        count = count + 1
    ENDIF

ENDFOR

newrarr =   newrarr[0:p-1]
newvrarr = newvrarr[0:p-1]
newnrarr = newnrarr[0:q-1]
newnvrarr = newnvrarr[0:q-1]
radarr = radarr[0:count-1]
colorarr = colorarr[0:count-1]

device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

;device, filename = '/n/Godiva1/jkrick/A3888/cmd/cmdnocolor.ps', /portrait, $
;  BITS=8, scale_factor=0.9 , /color
device, filename = '/n/Godiva1/jkrick/A3888/cmd/junk.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot,[0,0,0],[0,0,0],  thick = 3, psym = 2, xthick = 3, ythick = 3,charthick = 5,$
                        xrange = [15,27], yrange = [-1,1.5], ystyle = 1,xstyle = 1,$
	  ytitle = "V-r", $;title = 'cmd looking for cluster galaxies', $
                        xtitle = "r (mag)", charsize = 1.5
oplot,mgalaxy.rmag, mgalaxy.vr,  thick = 3, psym = 2,color = colors.darkgray
oplot, newrarr, newvrarr, psym = 5, color = colors.black, thick = 3
oplot, newnrarr,newnvrarr, psym = 6 ,color = colors.black, thick = 3

;oplot,r1, yfit1, thick = 3;, color = colors.yellow
;oplot, r1, yfit1 + sig1;, thick = 3, color = colors.yellow
;oplot, r1, yfit1 - sig1;, thick = 3, color = colors.yellow

xvals = findgen(11) +16
oplot, xvals, ((result(0))*xvals) + result(1), thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)+ sig1 , thick = 3;, color = colors.orange
oplot, xvals, ((result(0))*xvals) + result(1)- sig1 , thick = 3;, color = colors.orange


print, "sig1", sig1
;plothist,mgalaxy.vr,xhist,yhist,bin=0.05
print,"mean", mean(mgalaxy.vr)

;'----------------------------------------------------
newrad = sqrt((mgalaxy.xcen - 1556.)^2. + (mgalaxy.ycen - 1784.)^2.)
plot, newrad, mgalaxy.vr, thick = 3, psym = 2,yrange=[-1.0,1.5], color = colors.gray
oplot, radarr, colorarr, psym = 2, color = colors.red

coeff = ROBUST_LINEFIT( newrad, mgalaxy.vr, yfit, sig, coeff_sig)
err = dindgen(40) - dindgen(40) + 1
start = [-0.02,1.0]
r1j = fltarr(40)
y1j = fltarr(40)
r1j = newrad[3:42]
y1j = yfit[3:42]
sortindex = sort(r1j)
sortr1 = r1j[sortindex]
sortyfit1 = y1j[sortindex]
result = MPFITFUN('linear',sortr1, sortyfit1,err, start)
oplot, newrad, ((result(0))*newrad) + result(1), thick = 3, color = colors.black

coeff = ROBUST_LINEFIT( radarr, colorarr, yfit, sig, coeff_sig)
err = dindgen(40) - dindgen(40) + 1
start = [-0.02,1.0]
r1j = fltarr(40)
y1j = fltarr(40)
r1j = radarr[3:42]
y1j = yfit[3:42]
sortindex = sort(r1j)
sortr1 = r1j[sortindex]
sortyfit1 = y1j[sortindex]
result = MPFITFUN('linear',sortr1, sortyfit1,err, start)
oplot, radarr, ((result(0))*radarr) + result(1), thick = 3, color = colors.red

print, "cluster density", density, density/(!Pi*(1182^2))
print, "group density", densityg, densityg/(!Pi*(400^2))
print, "rsum, Vsum", rsum, Vsum
device, /close
set_plot, mydevice

undefine, rgalaxy
undefine, lgalaxy
undefine, mgalaxy
undefine, Vgalaxy

close, nolun
free_lun, nolun
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

