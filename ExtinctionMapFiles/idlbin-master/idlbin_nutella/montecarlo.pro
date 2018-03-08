pro montecarlo
!P.multi = [0,2,2]
numclusters = 7

;do the experiment 10 times
nexperiments = 1000

;now need the distance between each elements of randdist
randdist = fltarr(50000); fltarr(numclusters * (numclusters - 1))  ;there must be a smarter way to do this.
avgranddist = fltarr(nexperiments)
k = 0
count60 = 0
for l = 0, nexperiments -1 do begin
;   print, 'l new experiment', l
   randx = randomu(seed, numclusters)
   randy = randomu(seed, numclusters)
   randx = randx * 20.*60.
   randy = randy * 20.*60.
;   plot, findgen(1200), findgen(1200),/nodata
;   xyouts, randx, randy, 'x'
   for i = 0, numclusters-1 do begin
      for j = i + 1, numclusters -1 do begin
;         print, i, j, randx(i), randx(j), randy(i) , randy(j), sqrt((randx(i) - randx(j))^2 +  (randy(i) - randy(j))^2)
         randdist(k) = sqrt((randx(i) - randx(j))^2 +  (randy(i) - randy(j))^2)
         if randdist(k) le 108 then count60 = count60 + 1
         k = k + 1
      endfor
   endfor
;   print, "mean, sttdev randdist", mean(randdist[0:k-1]), stddev(randdist[0:k-1])
   
   avgranddist[l] = mean(randdist)
endfor

randdist = abs(randdist[0:k-1])

;plothist, randdist, bin = 1, xthick =1, yhtick = 1, charthick = 1, xrange=[0,2000]
print, 'mean(randdist), stddev(randdist)', mean(randdist), stddev(randdist)

nsigma = (60. - mean(randdist)) / stddev(randdist)
print, nsigma
PRINT, 'probability', (GAUSS_PDF( (60. - mean(randdist))/(stddev(randdist)) )  ) 

print, 'count 60, nexperiments', count60, nexperiments

;maybe I can ask how many times I get distances <=60
end
