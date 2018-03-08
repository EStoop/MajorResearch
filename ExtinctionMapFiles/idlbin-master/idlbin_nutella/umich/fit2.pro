
pro fit2


device, true=24
device, decomposed=0
colors = GetColor(/Load, Start=1)
close, /all		;close all files = tidiness


OPENR, lun,'/n/Sheriff1/jkrick/sep00/satlist', /GET_LUN
;read in the radial profile
rows= 9

unsat = FLTARR(rows)
sat = FLTARR(rows)

u = 0.0
s = 0.0

FOR j=0,rows-1 DO BEGIN
	READF, lun, u,s
	unsat(j) = u
	sat(j) = s
ENDFOR

close, lun
free_lun, lun

err = dindgen(rows) - dindgen(rows) + 1

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit

sortindex = Sort(unsat)
sortedunsat = unsat[sortindex]
sortedsat = sat[sortindex]


start = [0.9]

result = MPFITFUN('linear',sortedunsat, sortedsat, err, start, maxiter = 300)
;print, "sat", sat, "result*sat",result(0)*sat


;set up for plotting
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'


device, filename = '/n/Sheriff1/jkrick/sep00/satcounts.ps', /portrait, $
                BITS=8, scale_factor=0.9 , /color

plot, sortedunsat, sortedsat, psym = 2, thick=3, $
	XRANGE = [0,1.5E5], YRANGE = [0,1.5E5],$
	xtitle = 'counts per second from unsaturated stars',$	
	ytitle = 'counts per second from saturated stars', $
	title = 'is there a relation between saturated and unsaturated?', $
	subtitle ='blue = best fit line, red = slope of 1

oplot, sortedunsat, result(0)*sortedunsat, color = colors.blue

oplot, sortedunsat, sortedunsat, color = colors.red, linestyle = 1
device, /close
set_plot, mydevice
END
