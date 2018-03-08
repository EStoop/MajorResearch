PRO pixinten
close, /all
fits_read, "/n/Godiva1/jkrick/A3888/centerV2000.block.fits", data, header
;fits_read, "/n/Godiva1/jkrick/A3888/centersmo.fits", data, header
data = data - 2.0

;nominal centers by eye
xcenter = 1026.
ycenter = 1018.

;counters = intarr(14)    ;initialized to zero
ta = 0
a =0 
b = 0
c = 0
d = 0
e = 0
f = 0
g = 0
h = 0
i = 0
j = 0
k = 0
l = 0
m = 0
n = 0
o= 0
p = 0
q = 0
s = 0
t = 0
u = 0
v = 0
w = 0
z = 0
aa =0 
ab = 0
ac = 0
ad = 0
ae = 0
af = 0
ag = 0
ah = 0
ai = 0
aar = fltarr(5000)
bar = fltarr(5000)
car = fltarr(5000)
dar = fltarr(5000)
ear = fltarr(5000)
far = fltarr(5000)
gar = fltarr(5000)
har = fltarr(5000)
iar = fltarr(5000)
jar = fltarr(5000)
kar = fltarr(5000)
lar = fltarr(5000)
mar = fltarr(5000)
nar = fltarr(5000)
oar = fltarr(5000)
par = fltarr(5000)
qar = fltarr(5000)
sar = fltarr(5000)
tar = fltarr(5000)
uar = fltarr(5000)
var = fltarr(5000)
war = fltarr(5000)
zar = fltarr(5000)
aaar = fltarr(5000)
abar = fltarr(5000)
acar = fltarr(5000)
adar = fltarr(5000)
aear = fltarr(5000)
afar = fltarr(5000)
agar = fltarr(5000)
ahar = fltarr(5000)
aiar = fltarr(5000)
test = fltarr(5000)
;make annuli with similar number of pixels in each
FOR x = 0., 2000.,1 DO BEGIN
    FOR y = 0., 2000. , 1 DO BEGIN
        IF data[x,y] GT - 1 THEN BEGIN
            r = sqrt((abs(x - xcenter)^2) + (abs(y - ycenter)^2))
             IF r GE 70. AND r LT 100. THEN BEGIN
                test[ta] = data[x,y]
                ta = ta + 1
            ENDIF
             IF r GE 100. AND r LT 118. THEN BEGIN
                aar[a] = data[x,y]
                a = a + 1
            ENDIF
            IF r GE 118. AND r LT 131. THEN BEGIN
                bar[b] = data[x,y]
                b = b + 1
            ENDIF
            IF r GE 131. AND r LT 141. THEN BEGIN
                car[c] = data[x,y]
                c = c + 1
            ENDIF
            IF r GE 141. AND r LT 149. THEN BEGIN
                dar[d] = data[x,y]
                d = d + 1
            ENDIF
            IF r GE 149. AND r LT 156. THEN BEGIN
                ear[e] = data[x,y]
                e = e + 1
            ENDIF
            IF r GE 156. AND r LT 162. THEN BEGIN
                far[f] = data[x,y]
                f = f + 1
            ENDIF
            IF r GE 162. AND r LT 167. THEN BEGIN
                gar[g] = data[x,y]
                g = g + 1
            ENDIF
            IF r GE 167. AND r LT 171. THEN BEGIN
                har[h] = data[x,y]
                h = h + 1
            ENDIF
            IF r GE 171. AND r LT 175. THEN BEGIN
                iar[i] = data[x,y]
                i = i + 1
            ENDIF
            IF r GE 175. AND r LT 179. THEN BEGIN
                jar[j] = data[x,y]
                j = j + 1
            ENDIF
            IF r GE 179. AND r LT 183. THEN BEGIN
                kar[k] = data[x,y]
                k = k + 1
            ENDIF
            IF r GE 183. AND r LT 187. THEN BEGIN
                lar[l] = data[x,y]
                l = l + 1
            ENDIF
            IF r GE 187. AND r LT 191. THEN BEGIN
                mar[m] = data[x,y]
                m = m + 1
            ENDIF
            IF r GE 191. AND r LT 195. THEN BEGIN
                nar[n] = data[x,y]
                n = n + 1
            ENDIF
            IF r GE 195. AND r LT 199. THEN BEGIN
                oar[o] = data[x,y]
                o = o + 1
            ENDIF
            IF r GE 199. AND r LT 203. THEN BEGIN
                par[p] = data[x,y]
                p = p + 1
            ENDIF
            IF r GE 203. AND r LT 207. THEN BEGIN
                qar[q] = data[x,y]
                q = q + 1
            ENDIF
            IF r GE 207. AND r LT 211. THEN BEGIN
                sar[s] = data[x,y]
                s = s + 1
            ENDIF
            IF r GE 211. AND r LT 214. THEN BEGIN
                tar[t] = data[x,y]
                t = t + 1
            ENDIF
            IF r GE 214. AND r LT 217. THEN BEGIN
                uar[u] = data[x,y]
                u = u + 1
            ENDIF
            IF r GE 217. AND r LT 219. THEN BEGIN
                var[v] = data[x,y]
                v = v + 1
            ENDIF
            IF r GE 219. AND r LT 221. THEN BEGIN
                war[w] = data[x,y]
                w = w + 1
            ENDIF
            IF r GE 221. AND r LT 223. THEN BEGIN
                zar[z] = data[x,y]
                z = z + 1
            ENDIF
             IF r GE 223. AND r LT 225. THEN BEGIN
                aaar[aa] = data[x,y]
                aa = aa + 1
            ENDIF
            IF r GE 225. AND r LT 227. THEN BEGIN
                abar[ab] = data[x,y]
                ab = ab + 1
            ENDIF
            IF r GE 227. AND r LT 229. THEN BEGIN
                acar[ac] = data[x,y]
                ac = ac + 1
            ENDIF
            IF r GE 229. AND r LT 231. THEN BEGIN
                adar[ad] = data[x,y]
                ad = ad + 1
            ENDIF
            IF r GE 231. AND r LT 233. THEN BEGIN
                aear[ae] = data[x,y]
                ae = ae + 1
            ENDIF
            IF r GE 233. AND r LT 235. THEN BEGIN
                afar[af] = data[x,y]
                af = af + 1
            ENDIF
            IF r GE 235. AND r LT 237. THEN BEGIN
                agar[ag] = data[x,y]
                ag = ag + 1
            ENDIF
            IF r GE 237. AND r LT 239. THEN BEGIN
                ahar[ah] = data[x,y]
                ah = ah + 1
            ENDIF
            IF r GE 239. AND r LT 241. THEN BEGIN
               aiar[ai] = data[x,y]
                ai = ai + 1
            ENDIF

        ENDIF
    ENDFOR
ENDFOR

print,ta, a, b, c, d, e, f,g, h,i,j,k,l,m,n,o,p,q,s,t,u,v,w,z,aa, ab, ac, ad, ae, af,ag, ah,ai

test = test[0:ta-1]
aar = aar[0:a-1]
bar = bar[0:b-1]
car = car[0:c-1]
dar = dar[0:d-1]
ear = ear[0:e-1]
far = far[0:f-1]
gar = gar[0:g-1]
har = har[0:h-1]
iar = iar[0:i-1]
jar = jar[0:j-1]
kar = kar[0:k-1]
lar = lar[0:l-1]
mar = mar[0:m-1]
nar = nar[0:n-1]
oar = oar[0:o-1]
par = par[0:p-1]
qar = qar[0:q-1]
sar = sar[0:s-1]
tar = tar[0:t-1]
uar = uar[0:u-1]
var = var[0:v-1]
war = war[0:w-1]
zar = zar[0:z-1]
aaar = aaar[0:aa-1]
abar = abar[0:ab-1]
acar = acar[0:ac-1]
adar = adar[0:ad-1]
aear = aear[0:ae-1]
afar = afar[0:af-1]
agar = agar[0:ag-1]
ahar = ahar[0:ah-1]
aiar = aiar[0:ai-1]




;ps_open, file = "/n/Godiva1/jkrick/A3888/pixinten.ps", /portrait, /color, xsize = 6, ysize = 6
device, true=24
device, decomposed=0

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/Load, Start=1)

device, filename = '/n/Godiva1/jkrick/A3888/pixinten.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


plothist, test, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
plot, xhist, yhist, thick = 3, xrange = [-0.01,0.04], yrange = [0,3.5], charthick = 3, xthick = 3, $
ythick = 3, xtitle = "Intensity (counts/s)", title = "V-band Pixel Intensity Histograms", ystyle = 5

plothist, aar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist + 0.1, thick = 3

plothist, bar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist +0.2, thick = 3
plothist, car, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.3, thick = 3
plothist, dar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.4, thick = 3
plothist, ear, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.5, thick = 3
plothist, far, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.6, thick = 3
plothist, gar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.7, thick = 3
plothist, har, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.8, thick = 3
plothist, iar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+0.9, thick = 3
plothist, jar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.0, thick = 3
plothist, kar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.1, thick = 3
plothist, lar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.2, thick = 3
plothist, mar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.3, thick = 3
plothist, nar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.4, thick = 3

plothist, oar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.5, thick = 3
plothist, par, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.6, thick = 3
plothist, qar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.7, thick = 3
plothist, sar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.8, thick = 3
plothist, tar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+1.9, thick = 3
plothist, uar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.0, thick = 3
plothist, var, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.1, thick = 3
plothist, war, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.2, thick = 3
plothist, zar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.3, thick = 3

plothist, aaar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist +2.4, thick = 3
plothist, abar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist +2.5, thick = 3
plothist, acar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.6, thick = 3
plothist, adar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.7, thick = 3
plothist, aear, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.8, thick = 3
plothist, afar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+2.9, thick = 3
plothist, agar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+3.0, thick = 3
plothist,ahar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3
oplot, xhist, yhist+3.1, thick = 3
plothist, aiar, xhist, yhist, bin = 0.001, thick = 3, /noplot, peak = 0.3


;oplot, [median(aar),median(bar),median(car),median(dar),median(ear),median(far),median(gar),median(har),median(iar),median(jar),median(kar),median(lar),median(mar),median(nar),median(oar),median(par),median(qar),median(sar),median(tar),median(uar),median(var),median(war),median(zar),median(aaar),median(abar),median(acar),median(adar),median(aear),median(afar),median(agar),median(ahar),median(aiar)], [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0, 3.1, 3.2, 3.3], color = colors.blue, thick = 3

oplot, [mean(test),mean(aar),mean(bar),mean(car),mean(dar),mean(ear),mean(far),mean(gar),mean(har),mean(iar),mean(jar),mean(kar),mean(lar),mean(mar),mean(nar),mean(oar),mean(par),mean(qar),mean(sar),mean(tar),mean(uar),mean(var),mean(war),mean(zar),mean(aaar),mean(abar),mean(acar),mean(adar),mean(aear),mean(afar),mean(agar),mean(ahar),mean(aiar)], [0.27, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0, 3.1, 3.2, 3.3, 3.4], color = colors.red, thick = 5



;plot, [0,0,0],[0,0,0], xrange = [-0.01, 0.04], yrange = [0,1], charthick = 3, xthick = 3, ythick = 3, xtitle = "Intensity (counts/s)", title = "Pixel Intensity Histograms"
;plothist, aar, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.3
;plot, xhist, yhist, thick = 3, xrange = [-0.02, 0.05], yrange = [0,1]
;plothist, bar, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.4
;plothist, car, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.5
;plothist, dar, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.6
;plothist, ear, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.7
;plothist, far, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.8
;plothist, gar, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 0.9
;plothist, har, xhist, yhist, bin = 0.001, thick = 3, /overplot, peak = 1.0


;ps_close, /noprint, /noid

device, /close
set_plot, mydevice




END

