pro find_dogs
close, /all
restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'
!P.multi = [0,1,1]


redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

;needs to be detected at 24, K, and r
;having trouble with catalog matching 24 and r so just try 24 with any
;r (or undetected r) and the same with K
basic = where(objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90 and objectnew.rmag gt 0 and objectnew.wirckmag gt 0 )

;but fix those undetections to the limit in r
a = where(objectnew[basic].rmag gt 90)
objectnew[basic[a]].rmag = 26.6

;but fix those undetections to the limit in K
a = where(objectnew[basic].wirckmag gt 90)
objectnew[basic[a]].wirckmag = 21.2


flux24 = abtojansky(objectnew[basic].mips24mag) ;in jy
flux24 = flux24*1E6  ; in microjy

rflux = abtojansky(objectnew[basic].rmag)
rflux = rflux*1E6


dog = where( flux24 / rflux gt 1000 and flux24 gt 100 and objectnew[basic].rmag - objectnew[basic].wirckmag gt 4.5)

print, n_elements(dog)


;plothyperz, basic[dog], '/Users/jkrick/nutella/spitzer/dogs.ps'

;make *.reg file for viewing objcts in ds9
openw, outlunred, '/Users/jkrick/nutella/spitzer/dogs.reg', /get_lun
printf, outlunred, 'fk5'
for i=0, n_elements(dog) -1 do  printf, outlunred, 'circle( ', objectnew[basic[dog[i]]].ra, objectnew[basic[dog[i]]].dec, ' 3")'
close, outlunred
free_lun, outlunred


end
;vetting by eye

;number, objectnew number, reason

;1, 1372, is right next to satstar, so r was undetected
;2, 1956, underneath of the diffraction spike in r band so undetected
;#3, 2219, r band not deblended correctly. could do photometry here.
;*4, 2556, questionable match
;5, 2665, bright galaxy center
;6,2886, bright star
;#7, 2896, mismatch, there is an r detection there, could do photometry here.
;#8, 2938, mismatch with 2362, there is an r detection there, could do photometry here.
;#9, 2980, extended source mismatched, 
;*10, 2993, questionable match right next to saturated star
;11, 3020, bright star
;#12, 3029, mismatch, is a bright galaxy there
;13, 3226, edge of frame
;*14, 3227, mismatch, near bright star, but potentially there.
;#15, 3255, mismatch on a spiral
;#16, 3290, same object as 3255
;17, 3293, bright galaxy
;18, 3298, same big bright galaxy
;#19, 3304, mismatch
;#20, 51064, mismatch
;#21, 56979, mismatch
;#22, 57234, mismatch
;#23, 57853,mismatch
;#24, 57984, mismatch
;25, 58645, mismatch
;26, 61378, mismatch
;27, 62148, mismatch
;28, 63712, mismatch, cool galaxy
;29, 76941, mismatch
;30, 77964, mismatch
;31, 86627, mismatch
