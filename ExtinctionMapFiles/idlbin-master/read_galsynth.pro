pro read_galsynth
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

;filename = '/Users/jkrick/Downloads/job_1000001571/result.0'

;-----------------------------------


filename = '/Users/jkrick/Virgo/grasil/elli1_01.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

plot, lambda, total * lambda, xrange = [ 1E3,1E8] ,/xlog, /ylog, yrange=[1E10,1E18]
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3

;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_02.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = redcolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3


;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_04.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = bluecolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3

;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_08.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = greencolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3

;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_1_5.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = yellowcolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3

;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_5.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = cyancolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3

;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_8.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = orangecolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3

;-----------------------------------
filename = '/Users/jkrick/Virgo/grasil/elli1_13.0.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"

oplot, lambda, total * lambda, color = purplecolor
;oplot,lambda,  cirrus* lambda, linestyle = 1
;oplot,lambda,  star1* lambda, linestyle = 2
;;;oplot,lambda,  star2* lambda, linestyle = 3
;oplot,lambda,  mc* lambda, linestyle = 4
;;;oplot,lambda,  star3* lambda, linestyle = 5
;;;oplot,lambda,  sub1* lambda, linestyle = 2
;;;oplot,lambda,  sub2* lambda, linestyle = 3
;-----------------------------------


end
