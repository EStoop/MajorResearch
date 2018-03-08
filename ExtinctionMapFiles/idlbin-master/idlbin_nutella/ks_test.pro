pro ks_test

ps_open, filename='/Users/jkrick/nep/clusters/kstest.ps', /portrait, /square, xsize=6, ysize=4, /color
!p.multi=[0,2,1]
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

;set up junk data to test
data1 = randomn(seed, 40)
data2 = randomn(seed, 70)

;data 1 is redshift 1.0 clusters
restore, '/Users/jkrick/idlbin/object.old.sav'
restore, '/Users/jkrick/nep/clusters/memberarr.txt'
dm1 = 41.1  ;including all stuff already
data1 = objectnew[memberarr].irac1mag - dm1
;only look at the distributions betwen K = -18 and -22
a = where(data1 lt -18 and data1 gt -22)
data1 = data1(a)

;data 2 is redshift 0.1 clusters
restore, '/Users/jkrick/nep/clusters/cosmos_allmemberk.txt'
dm2 = 38.2
data2 = allmemberk - dm2
b = where(data2 lt -18 and data2 gt -22)
data2 = data2(b)


kstwo, data1, data2, D, prob
print, 'd, prob', d, prob

;want to plot the two distributions
sortdata1 = data1[sort(data1)]
sortdata2 = data2[sort(data2)]

N1 = n_elements(sortdata1)
N2 = n_elements(sortdata2)

f1 = (findgen(N1) + 1.)/ N1
f2 = (findgen(N2) + 1.)/ N2

plot, sortdata1, f1, psym = 10, xtitle='Mk', ytitle='cumulative distribution', thick=3, xthick=3, ythick=3, charthick=3
oplot, sortdata2, f2, psym = 10, linestyle = 1, thick=3
xyouts, -21, 0.8, 'clusters'



;------------------------------------------------
;now consider what the background distributions look like.

;backdata 1 is IRAC dark field background
restore, '/Users/jkrick/nep/clusters/randmemberarr.txt'
backdata1 = objectnew[randmemberarr].irac1mag - dm1
c = where(backdata1 lt -18 and backdata1 gt -22)
backdata1 = backdata1(c)


;backdata2 is COSMOS background
restore, '/Users/jkrick/nep/clusters/cosmos/backmemberarr.txt'
backdata2= comosbackmemberarr - dm2
d = where(backdata2 lt -18 and backdata2 gt -22)
backdata2 = backdata2(d)

kstwo, backdata1, backdata2, D, prob
print, 'background d, prob', d, prob

;want to plot the two distributions
sortbackdata1 = backdata1[sort(backdata1)]
sortbackdata2 = backdata2[sort(backdata2)]

N3 = n_elements(sortbackdata1)
N4 = n_elements(sortbackdata2)

f3 = (findgen(N3) + 1.)/ N3
f4 = (findgen(N4) + 1.)/ N4

plot, sortbackdata1, f3, xtitle='Mk', ytitle='cumulative distribution', psym = 10, thick=3, xthick=3, ythick=3, charthick=3
;, color=bluecolor
oplot, sortbackdata2, f4, psym = 10, linestyle = 1, thick=3;, color=bluecolor
xyouts, -21, 0.8, 'background'
ps_close, /noprint,/noid

end
