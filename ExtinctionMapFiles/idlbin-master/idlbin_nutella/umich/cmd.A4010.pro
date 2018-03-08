PRO cmd
close,/all
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; february 2004
;
;this code makes a color magnitude diagram of "galaxies" in A2556
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fits_read, "/n/Godiva4/jkrick/A4010/original/fullr.wcs.fits", data, header
;fits_read, "/n/Godiva3/jkrick/A2556/original/fullr.fits", data, header

numoobjects = 45000
rgalaxy = replicate({object, xcenter:0D,ycenter:0D,fluxa:0D,maga:0D,fluxb:0D,magb:0D,fwhm:0D,isoarea:0D},numoobjects)
Vgalaxy = replicate({Vobject, Vxcenter:0D,Vycenter:0D,Vfluxa:0D,Vmaga:0D,Vfluxb:0D,Vmagb:0D,Vfwhm:0D,Visoarea:0D},numoobjects)
mgalaxy = replicate({what, xcen:0D,ycen:0D, vr:0D,rmag:0D, member:0D},numoobjects)

i = 0
rsum = 0.
Vsum = 0.

openr, lun, "/n/Godiva4/jkrick/A4010/cmd/SExtractor.r.cat", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF,lun,junk, rx, ry, fluxaper,magaper,erraper,fluxbest,magbest,errbest,fwhm,rra,rdec,iso
    IF (fwhm GT 4.2) AND (fluxbest GT 0.0) AND (iso GT 6.0) THEN BEGIN   
      ;no bother including the stars and junk
        rgalaxy[i] = {object,  rx, ry, fluxaper,magaper,fluxbest,magbest,fwhm,iso}
        i = i + 1
    ENDIF
ENDWHILE

rgalaxy = rgalaxy[0:i - 1]
print, "there are i,",i," r galaxies"
close, lun
free_lun, lun

;read in membership information
numogals = 212
l = 0
lgalaxy = replicate({couch,x:0D,y:0D},numogals)
openr, lun3, "/n/Godiva4/jkrick/A4010/cmd/collins.4010.member.txt", /GET_LUN
WHILE (NOT eof(lun3)) DO BEGIN
    READF, lun3,junk,rh,rm,rs,dd,dm,ds,junk,junk,junk
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.);
    adxy, header, r, d, xcen, ycen
;    print, r, d, xcen,ycen
    lgalaxy[l] = {couch,xcen,ycen}
    l = l +1

ENDWHILE
close, lun3
free_lun, lun3
lgalaxy[l] = {couch, 1513.76,1729.38  }
lgalaxy = lgalaxy[0:l]

;print, lgalaxy

;read in membership information
numogals = 212
n = 0
mn = 0
l2galaxy = replicate({couch,x:0D,y:0D},numogals)
l3galaxy = replicate({couch,x:0D,y:0D},numogals)
openr, lun3, "/n/Godiva4/jkrick/A4010/cmd/katgert.member.txt", /GET_LUN
WHILE (NOT eof(lun3)) DO BEGIN
    READF, lun3,junk,rh,rm,rs,dd,dm,ds,junk,cz,junk
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
;    print, r, d, xcen,ycen
    IF (cz GT 2.7E4) AND (cz LT 3.05E4) THEN BEGIN  ;member
        l2galaxy[n] = {couch,xcen,ycen}
        n = n + 1
    ENDIF ELSE BEGIN  ;nonmember
        l3galaxy[mn] = {couch,xcen,ycen}
        mn = mn + 1
    ENDELSE
ENDWHILE
close, lun3
free_lun, lun3
l2galaxy = l2galaxy[0:n-1]
l3galaxy = l3galaxy[0:mn-1]

;print, "l3galaxy",l3galaxy
;read in the V magnitudes and match them up to make a V-r

m = 0
ship = 0
nship = 0
vra = 0D
vdec = 0D
xshift = 0.;-190
yshift = 0.;+60
openr, lun2, "/n/Godiva4/jkrick/A4010/cmd/SExtractor.B.cat", /GET_LUN

WHILE (NOT eof(lun2)) DO BEGIN
    READF, lun2,junk, vx, vy, vfluxaper,vmagaper,verraper,vfluxbest,vmagbest,$
     verrbest,vfwhm, vra,vdec,viso
;    xyad,header,vx,vy,vra,vdec     ;find the ra and dec from the better astrometry image
    membership = 0.                  ;assume it is not a member

    IF (vfwhm GT 3.9) AND (vfluxbest GT 0.0) AND (viso GT 6.0)THEN BEGIN ;not star,junk
        Vgalaxy[i] = {Vobject,  vx, vy, vfluxaper,vmagaper,vfluxbest,vmagbest,vfwhm,viso}
        FOR j = 0, i - 1, 1 DO BEGIN

            IF j LT l  THEN BEGIN      ;is this v galaxy really a member          
                IF (vx LT lgalaxy[j].x + 3.0) AND (vx GT lgalaxy[j].x - 3.0) AND $
                     (vy LT lgalaxy[j].y + 3.0) AND (vy GT lgalaxy[j].y - 3.0) THEN BEGIN
                    membership =  1.
;;                    ship = ship + 1
                ENDIF
            ENDIF
            IF j LT n  THEN BEGIN      ;is this v galaxy really a member          
                IF (vx LT l2galaxy[j].x + 5.0) AND (vx GT l2galaxy[j].x - 5.0) AND $
                     (vy LT l2galaxy[j].y + 5.0) AND (vy GT l2galaxy[j].y - 5.0) THEN BEGIN
                    membership =  1.
;                    ship = ship + 1
                ENDIF

            ENDIF
            IF j LT mn  THEN BEGIN      ;is this v galaxy really a nonmember          
                IF (vx LT l3galaxy[j].x + 5.0) AND (vx GT l3galaxy[j].x - 5.0) AND $
                     (vy LT l3galaxy[j].y + 5.0) AND (vy GT l3galaxy[j].y - 5.0) THEN BEGIN
                    membership =  10.
;                    ship = ship + 1
                ENDIF

            ENDIF

            IF (vx LT rgalaxy[j].xcenter + 1.5 +xshift) AND $
              (vx GT rgalaxy[j].xcenter - 1.5 + xshift) AND $
              (vy LT rgalaxy[j].ycenter + 1.5 + yshift) AND $ 
              (vy GT rgalaxy[j].ycenter - 1.5 + yshift)  THEN BEGIN
                
                mgalaxy[m] = {what, vx,vy,vmagbest - rgalaxy[j].magb ,rgalaxy[j].magb,membership}

;                IF rgalaxy[j].magb LT 16.0 AND $
;                  (vmagbest - rgalaxy[j].magb GT 0.0)  THEN  $
;                  print, "brightblue", vx,vy, rgalaxy[j].magb, $
;                     vmagbest - rgalaxy[j].magb

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
    IF (mgalaxy[se].rmag LE 20.) AND (mgalaxy[se].rmag GT 16.0) $
      AND (mgalaxy[se].vr LT 2.3) AND (mgalaxy[se].vr GT 1.4) THEN BEGIN
        r1(t) = mgalaxy[se].rmag
        vr1(t) = mgalaxy[se].vr
        t = t +1
    ENDIF 
    IF (mgalaxy[se].rmag LT 22.) AND (mgalaxy[se].rmag GT 16.0) THEN BEGIN
        r2(u) = mgalaxy[se].rmag
        vr2(u) = mgalaxy[se].vr
        u = u + 1
    ENDIF

ENDFOR
print, t, " the number of galaxies < 21"
r1 = r1[0:t-1]
vr1 = vr1[0:t-1]
r2 = r2[0:u-1]
vr2 = vr2[0:u-1]

;biweight fit (robust for non-gaussian distributions)
;coeff = ROBUST_LINEFIT( mgalaxy.rmag, mgalaxy.vr, yfit, sig, coeff_sig)
coeff1 = ROBUST_LINEFIT(r1,vr1, yfit1, sig1, coeff_sig1)
coeff2 = ROBUST_LINEFIT(r2,vr2, yfit2, sig2, coeff_sig2)

print, "sig1, sig2", sig1,sig2
err = dindgen(40) - dindgen(40) + 1
start = [-0.02,1.0]
r1j = fltarr(40)
y1j = fltarr(40)
r1j = r1[3:42]
y1j = yfit1[3:42]
sortindex = sort(r1)
sortr1 = r1[sortindex]
sortyfit1 = yfit1[sortindex]
result = MPFITFUN('linear',sortr1, sortyfit1,err, start)
print, "result", result

;--------------------------------------------------------
;for every galaxy , given it's r value, does it's v-r value lie
;within the acceptable range?
;-----------------------------------------------------------
openw, nolun, "/n/Godiva4/jkrick/A4010/cmd/member.txt",/get_lun
memgalincut = 0.
nmemgalincut = 0.
memgalnincut = 0.
nmemgalnincut = 0.
galcount = 0.
num = 0
density=0
memarray= fltarr(20000)
memcount = 0
xarr = fltarr(20000)
yarr = fltarr(20000)


FOR a = 0, m-1 ,1 DO BEGIN

    IF (mgalaxy[a].member GT 0 ) AND ( mgalaxy[a].member LT 9.)THEN num = num + 1

     ;is v-r LT yfit[a] + sig and GT yfit[a] - sig
    limithigh = ((result(0))*mgalaxy[a].rmag) + result(1)+ sig1 
    limitlow = ((result(0))*mgalaxy[a].rmag) + result(1)- sig1 
    IF (mgalaxy[a].vr LT limithigh )AND (mgalaxy[a].vr GT limitlow) THEN BEGIN
       ;galaxy "a" makes the color cut
        galcount = galcount + 1
        printf, nolun, mgalaxy[a].xcen ,mgalaxy[a].ycen , mgalaxy[a].rmag, $
          mgalaxy[a].vr + mgalaxy[a].rmag
        memarray[memcount] = mgalaxy[a].rmag
        xarr[memcount] = mgalaxy[a].xcen
        yarr[memcount] = mgalaxy[a].ycen

        memcount = memcount + 1


        IF (mgalaxy[a].member GT 0.) AND ( mgalaxy[a].member LT 9.)THEN BEGIN
            memgalincut = memgalincut + 1
        ENDIF

        IF (mgalaxy[a].member GT 9.) THEN BEGIN
            nmemgalincut = nmemgalincut + 1
        ENDIF

        ;count the flux in member galaxies

        dist1 = sqrt((mgalaxy[a].xcen-1514.)^2 + (mgalaxy[a].ycen - 1729.)^2) 
;        dist2 = sqrt((mgalaxy[a].xcen-1157.)^2 + (mgalaxy[a].ycen - 2267.)^2) 

        IF dist1 LE 997. THEN BEGIN
              rsum = rsum + 10^((22.04 - mgalaxy[a].rmag)/2.5)
              Vsum = Vsum + 10^((22.19 - (mgalaxy[a].vr + mgalaxy[a].rmag))/2.5)
          ENDIF
          IF dist1 LE 1034 AND mgalaxy[a].rmag LT 19.67 THEN density = density + 1
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
;#######################################
print, "percent of members included in color cut", memgalincut/(memgalincut+memgalnincut) *100.
;print, "percent of nonmembers included in color cut", nmemgalincut/(nmemgalnincut + nmemgalincut) *100.
print, "percent of all objects that made the color cut", galcount / m *100.

print,"memgalincut", memgalincut
print,"memgalnincut", memgalnincut
;print,"nmemgalincut", nmemgalincut
;print,"nmemgalnincut", nmemgalnincut
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
;    dist1 = sqrt((xarr(num)-1514.)^2 + (yarr(num) - 1729.)^2) 
;    IF (memarray(num) LT ( memarray[9] + 3)) AND (dist1 LE 1112.) THEN density = density + 1
;ENDFOR
;-------------------------------------------------------------
;for plotting up the confirmed members and non-members
;----------------------------------------------------------------
newrarr = fltarr(200)
newvrarr = fltarr(200)
newnrarr = fltarr(200)
newnvrarr = fltarr(200)
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
device, filename = '/n/Godiva4/jkrick/A4010/cmd/cmd.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot,[0,0,0],[0,0,0],  thick = 3, psym = 2, xthick = 3, ythick = 3,charthick = 5,$
          xrange = [14,25], yrange=[0,3.5],ystyle = 1,xstyle = 1,$
	  ytitle = "B-r", $;title = 'cmd looking for cluster galaxies', $
                        xtitle = "r (mag)", charsize = 1.5
oplot,mgalaxy.rmag, mgalaxy.vr,  thick = 3, psym = 2,color = colors.darkgray
oplot, newrarr, newvrarr, psym = 5, color = colors.black, thick = 3
oplot, [15.116900],[1.8406000],psym = 5, color = colors.black, thick = 3

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
device, /close
set_plot, mydevice

undefine, rgalaxy
undefine, lgalaxy
undefine, mgalaxy
undefine, Vgalaxy


print, "rsum,Bsum ", rsum, Vsum
print, "density", density
close, nolun
free_lun, nolun
END
