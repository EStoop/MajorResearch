pro virgo_find_cluster_2
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
;ps_start, filename= '/Users/jkrick/virgo/redshift_dist.ps'

;a large region including the potential cluster
filename = '/Users/jkrick/virgo/virgo_sdss.csv'
readcol, filename, ra, dec, u, err_u, g, err_g, r ,err_r, i ,err_i,z ,err_z, photoz, photozerr, photozcc2, photozerrcc2,photozd1, photozerrd1

a = where( photozerrcc2 lt 0.2)
plothist, photozcc2(a), bin = 0.05, thick = 3,/noplot

;plot the cluster histogram
clusterra = 186.76731
clusterdec = 13.146311

good = where(sqrt((ra - clusterra)^2 + (dec-clusterdec)^2) lt .04 and photozerrcc2 lt 0.2,ngood  )
print, 'there are', ngood, ' fitting objects'
clusterz = photozcc2(good)
;print, clusterz
plothist, clusterz, xhist, yhist,bin = 0.05,/halfbin, xrange = [0, 1.0], xtitle='SDSS Photometric Redshift', $
          ytitle = 'Number of Galaxies', yrange =[0,14], thick = 3, xthick = 3, ythick = 3, charthick = 3
 
;-------------------------------------
; what is the distribution of some random regions
;also need to know the standard deviation in each redshift bin

randra=[186.90415,186.78176,186.67262,186.91545,187.04347,187.10362,186.78735,186.87598,187.04522]
randdec=[13.2659, 13.236533,13.15763,13.17242,13.179726,13.05505,13.033075,12.864473,12.919441]

randz = fltarr(9,21)
;r = 0
for cand=0, n_elements(randdec) - 1 do begin
   good = where(sqrt((ra - randra(cand))^2 + (dec-randdec(cand))^2) lt .04 and photozerrcc2 lt 0.2,ngood  )
   plothist, photozcc2(good), xhist, yhist, bin = 0.05, /noplot, /halfbin
   yhist = histogram(photozcc2(good), bin =0.05, min = 0, max = 1)
;   print, 'yhist', yhist
   for i = 0, n_elements(yhist) -1 do randz(cand,i) = yhist(i)
;   print, 'randz', randz(cand,*)
endfor 
;print, 'test 0', randz[0,*]

;now find an average and a stddev
x = findgen(21) * 0.05
y = x
e = x
for i = 0, n_elements(x) - 1 do begin
   y(i) = mean(randz(*,i))
   e(i) = stddev(randz(*,i))
endfor

oploterror, x,y,e, psym = 10, linestyle = 2, thick = 3, errthick = 3
print, 'x', x
print, 'e', e

;is this cluster already detected in the sdss cluster catalogs?
readcol, '/Users/jkrick/virgo/sdss_koester_cat.txt', RAJ2000,	DEJ2000,	RAJ2000,	DEJ2000,	zph	,zsp,	LBr,	LBi,	LTr,	LTi,	Ngal,	NR200

sortkey = sort(raj2000)
sortra = raj2000(sortkey)
sortdec = dej2000(sortkey)
sortzph = zph(sortkey)
sortzsp = zsp(sortkey)

a = where(sortra gt 186.5 and sortra lt 187 and sortdec gt 13 and sortdec lt 13.3,ngood)
print, ngood
for i = 0, ngood - 1 do print, sortra(ngood(i)), sortdec(ngood(i)), sortzph(ngood(i)), sortzsp(ngood(i))

;ps_end
end
