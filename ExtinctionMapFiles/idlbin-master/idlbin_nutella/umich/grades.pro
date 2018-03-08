PRO grades

close, /all
;for sections 001 & 002
grade2 = fltarr(21)
grade2 = [87.89,77.16,86.61,87.45,67.73,92.32,89.64,82.62,92.16,95.58,93.68,71.35,91.16,92.54,88.88,86.93,82.76,76.14,72.14,86.04,72.95]

final2 = fltarr(21)
final2 = [0.89,0.81,0.69,0.74,0.81,0.93,0.86,0.84,0.95,0.97,0.91,0.54,0.86,0.96,0.81,0.85,0.75,0.72,0.7,0.88,0.66]


grade1 = fltarr(22)
grade1 = [76.06,84.66,49.41,83.38,90.74,93.88,85.94,94.09,89.25,75.86,94.14,84.32,60.54,66.37,91.76,85.99,62.67,88.24,92.94,73.11,86.53,88.08]

final1 = fltarr(22)
final1 = [0.76,0.86,0.57,0.68,0.89,0.97,0.89,0.99,0.81,0.71,0.92,0.83,0.73,0.8,0.99,0.82,0.63,0.92,0.95,0.79,0.95,0.89]


;for sections 003 & 004

grade3 = fltarr(22)
grade3 = [92.52,84.44,88.46,93.55,95.97,91.54,91.6,93.81,86.22,96.07,87.43,91.47,81.3,90.73,95.12,95.97,67.66,84.62,95.3,79.53,95.28,83.82]

final3 = fltarr(22)
final3 = [0.83,0.93,0.87,0.96,0.91,0.88,0.85,0.97,0.89,0.96,0.82,0.88,0.76,0.79,0.93,1,0.83,0.89,0.93,0.68,0.96,0.86]

grade4 = fltarr(22)
grade4 = [86.03,81.03,86.6,75.66,89.21,97.02,93.51,79.69,93.36,80.77,90.85,77.43,73.84,91.6,90.02,81.93,94.62,88.29,91.29,87.93,94.93,96.8]

final4 = fltarr(22)
final4 = [0.79,0.88,0.69,0.95,0.95,0.99,0.94,0.97,0.94,0.95,0.96,0.7,0.82,0.95,0.91,0.79,0.94,0.89,0.84,0.7,0.89,0.93]

;---------------
;f05
;-------------
final1 = fltarr(18)
final1 = [94.64,85.23,89.11,90.98,99.58,89.76,92.78,95.84,96.88,88.3,92.47,79.28,87.02,89.44,94.47,90.55,96.52,86.05]


final2 = fltarr(17)
final2 = [88.54,92.26,90.61,90.,90.36,90.4,85.03,95.68,93.54,93.67,94.22,96.09,95.,91.37,91.5,92.93,84.97]


final3 = fltarr(22)
final3 = [67.45,73.38,73.08,75.12,93.12,89.91,82.04,86.54,78.29,43.06,80.12,78.09,28.44,88.59,79.71,74.2,85.23,82.59,77.8,77.46,79.68,76.57]

final4 = fltarr(25)
final4 = [85.06,79.84,89.23,87.11,91.88,86.3,73.48,89.65,93.2,72.23,91.8,90.54,95.8,95.5,69.27,79.53,87.7,89.19,79.85,76.98,83.09,89.56,94.39,87.36,83.87]


overallfinal = [final3, final4]
print, mean(overallfinal), "  mean"
print, stddev(overallfinal), " stddev"
print, median(overallfinal), " median"
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva7/jkrick/teaching/127.2f05.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plothist,overallfinal, bin = .1, thick = 3,xthick=3,ythick=3, chartick = 3,xrange=[80,100]

device, /close
set_plot, mydevice

END
