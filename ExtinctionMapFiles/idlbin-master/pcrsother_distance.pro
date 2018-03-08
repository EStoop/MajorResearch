pro pcrsother_distance
;what is the distance between the catalog peakup star and the actual target

;for pcrsother test
;BD+68_968
;pcrsra = 17 57 29.17
;pcrsdec = 68 01 15.1
pcrsra = 269.37154167      
pcrsdec = 68.02086111

;BD +67_1044
;targetra = 17 58 54.69
;targetdec = 67 47 36.9
targetra = 269.72787500      
targetdec = 67.79358333


d = sphdist(pcrsra, pcrsdec, targetra, targetdec,/degrees)
print, d


;------------------------------
;for snapshot tests
;------------------------------
;pcrsra = 17 26 41.69
;pcrsdec = 59 44 55.5
pcrsra= 261.67370833
pcrsdec = 59.74875000

;HD158460
;targetra =17 25 41.354
;targetdec = 60 02 54.23
targetra =261.42229167      
targetdec =60.04838889

d = sphdist(pcrsra, pcrsdec, targetra, targetdec,/degrees)
print, d


end
