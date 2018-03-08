PRO correlations
close,/all

device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness


;mydevice = !D.NAME
;!p.multi = [0,2, 2]
!p.multi =[0,0,1]
;SET_PLOT, 'ps'

;device, filename = '/n/Godiva2/jkrick/correlations.ps', /portrait, $
;  BITS=8, scale_factor=0.9 , /color
ps_open, file = "/Users/jkrick/umich/icl/correlations.2.ps", /portrait, /color, ysize = 5

;cluster properties
ngals=[76,59,89,100,95,206,158,193,125,243]
rvir=[.86,.82,.81,.87,1.03,1.22,1.16,1.25,1.17,1.14]
density = ngals / (!PI*rvir^2)
;density2 = [75,56,88,89,72,142,92,103,104,150] ;dl
density2=[76,62,99,104,93,189,151,185,220,288]
mass=[2.82E14,8.3E14,2.49E14,25E14,2.7E14,25.5E14,31E14,18.9E14,26.3E14,38E14]
masserr=[0.37E14,1.93E14,0.89E14,11E14,1.1E14,10.8E14,10E14,11.1E14,1E14,20E14]
cd = [1,1,1,1,1,0,1,0,1,0]
rich = [1,0,1,1,1,2,2,3,2,3]
z = [0.048,0.058,0.062,0.087,0.096,0.15,0.18,0.23,0.31,0.31]
sigma=[845,827,628,1247,625,1102,0,1130,1388,1947]
sigmaerr=[280,120,61,249,127,137,0,150,128,292]
m1m2=[0.91,0.52,0.51,0.93,0.7,0.12,0.57,0.42,0.28,0.23]
m1m3=[1.05,0.55,0.62,1.11,0.72,0.17,0.64,0.56,0.53,0.24]

;k corrections, evolutionary corrections
;and a color correction for V -> B  all in magnitudes
rzke = [-0.002,-0.007,-0.009,-0.018,-0.018,-0.017,-0.001,-0.020,0.048,0.048]
Bzke = [0.189,0.218,0.233,0.361,0.369,1.456,1.584,1.753,1.994,1.994]

;galaxies flux within 1/3 rvir as measured by the CMD
;with bolometric ab.mag. of sun
;clusterrflux=[1.52E12,1.09E12,1.55E12,1.66E12,1.4E12,3.3E12,2.1E12,2.8E12,1.7E12,4.9E12]
;clusterBflux1=[2.97E11,2.8E11,2.53E11,2.8E11,2.6E11,2.6E12,1.7E12,2.02E12,1.13E12,2.6E12]

;with band specific ab.mag. of sun
;;dl 0.3 rvir
;clusterrflux=[1.37316e12,9.88316e11,1.39871e12,1.50285e12,1.26340e12,2.98057e12,1.88819e12,2.49043e12,1.53444e12,4.43162e12]

;clusterBflux=[5.81903e11,5.47446e11,4.95660e11,5.52095e11,5.10032e11,2.80083e12,1.78107e12,2.17948e12,1.21830e12,2.79093e12]

;da 0.25rvir
clusterBflux=[4.98576e+11,4.61669e+11,4.17017e+11,4.62209e+11,4.90112e+11,2.77008e+12,1.88308e+12,2.71212e+12,1.46880e+12,3.38012e+12]

clusterrflux=[1.17840e+12,8.51462e+11,1.18061e+12,1.24525e+12,1.21682e+12,2.94861e+12,1.99343e+12,3.10474e+12,1.84032e+12,4.61305e+12]

;after BO correction

;clusterBflux=[1.45*4.98576e+11,1.5*4.61669e+11,1.51*4.17017e+11,1.61*4.62209e+11,1.65*4.90112e+11,1.9*2.77008e+12,2.0*1.88308e+12,2.2*2.71212e+12,2.5*1.46880e+12,2.5*3.38012e+12]

;clusterrflux=[1.45*1.17840e+12,1.5*8.51462e+11,1.51*1.18061e+12,1.61*1.24525e+12,1.65*1.21682e+12,1.9*2.94861e+12,2.0*1.99343e+12,2.2*3.10474e+12,2.5*1.84032e+12,2.5*4.61305e+12]


;galaxies flux corrected for k correction, and ev. correction
clusterrzke=clusterrflux* 10^(-0.4*rzke)
clusterBzke=clusterBflux* 10^(-0.4*Bzke)


clusterrzkeerr=0.3*clusterrzke
clusterBzkeerr=0.3*clusterBzke
print, "clusterrzkerr", clusterrzkeerr,"clusterBzkeerr",clusterBzkeerr

;luminosity of the ICL as measured with circ. integration at the z of
;the cluster
lumr =[2.45E11,1.9E11,1.99E11,1.36E11,6.4E11,4.13E11,4.9E11,3.8E11,4.36E11,8.1E11]
lumrerr=[1.18E11,6.5E10,0.335E11,1.18E11,1.41E11,1.95E11,2.3E11,9.5E10,7.8E10,1.13E11]
lumrerrfactor=lumrerr/lumr

;lumB =[4.84E10,1.5E10,2.21E10,1.83E10,1.1E11,3.15E11,5.1E11,1.6E11,3.49E11,3.9E11]
lumBerrfactor=[0.21,0.53,0.68,1.05,0.36,.29,.34,.32,.21,.26]
;lumBerr = lumB * lumBerrfactor

;print, "lumrerrfactor,lumBerrfactor", lumrerrfactor, lumBerrfactor

;total ICL luminosities with elliptical profiles done manually
;from totalprof.pro

;with bolometric ab.mag. of sun
;lumr=[8.1798513e+10,1.3263833e+11,1.1856016e+11,5.7033328e+10,2.4979011e+11,4.0032872e+11,1.9496604e+11,3.8E11,1.2250008e+11,8.1E11]
;lumB=[2.2755918e+10,8.4698513e+09,8.7994546e+09,4.4883498e+09,4.1140539e+10,2.2649573e+11,2.2265260e+11,1.6E11,9.5292021e+10,3.9E11]

;with band specific ab. mag. of sun

; summing to the mask edge
;lumr = [7.3917206e+10,1.1985859e11,1.0713686e11,5.1538153e10,2.2572277e11,3.6175695e11,1.7618101e11,2.69e11  ,1.4633013e11,6.32e11 ]

;lumB = [4.4575311e+10  ,1.6591123e10  ,1.7236765e10  ,8.7919803e09
;,8.0587926e10  ,2.4381487e11  ,2.3967787e11  ,1.25e11
;,1.7410223e11  ,3.62e11] 

;summing all the way into the center

lumB =[  1.3625637e+11 ,  5.4120878e+10,  8.6240063e+10,  1.9132582e+10,  1.0836237e+11,  3.2734241e+11,  2.6760109e+11,1.71438e+11    , 2.3684290e+11 ,4.22928e+11    ]

lumr =[3.4204715e+11,1.3436338e+11,2.7464027e+11,7.4731070e+10,3.1468707e+11,4.3235596e+11,2.1670132e+11,3.44589e+11  ,2.3441482e+11,7.29485e+11  ]


;lumr =[3.4204715e+11 + 0.43*3.4204715e+11,1.3436338e+11 + 0.43*1.3436338e+11,2.7464027e+11 +0.43*2.7464027e+11,7.4731070e+10 + 0.43*7.4731070e+10,3.1468707e+11 + 0.43*3.1468707e+11,4.3235596e+11,2.1670132e+11+0.43* 2.1670132e+11,3.44589e+11  ,2.3441482e+11 +0.43*2.3441482e+11,7.29485e+11  ]

;lumr = lumr*rellfactor
;lumB = lumB*Bellfactor


;total ICL luminosities with elliptical profiles, corrected for flux
;dimming, k corrections, ev. corrections, and color for the 4 high z's
;ie total ICL flux for all clusters at z=0, and their errors
lumrzke = lumr * 10^(-0.4*rzke)
lumrzkeerr=lumrerrfactor*lumrzke

lumBzke=lumB * 10^(-0.4*Bzke)
lumBzkeerr = lumBerrfactor*lumBzke


print, "lumr", lumr
print, "lumB", lumB
print, "lumrzke", lumrzke
print, "lumBzke", lumBzke
print, "lumrzkeer,", lumrzkeerr, "lumBzkeerr", lumBzkeerr

;old fractions before elliptical profiles
iclfracV=[13.,5.,8.,8.,30.,13.,20.,6.,17.,9.]
iclfracVerr=[5.,3.,6.,9.,14.,4.,9.,3.,4.,4.]
iclfracr=[13.,14.,10.,10.,31.,13.,16.,11.,16,11.]
iclfracrerr=[7.,6.,3.,9.,12.,5.,9.,4.,3.,4.]
rfracerr = iclfracrerr / iclfracr
Bfracerr = iclfracVerr / iclfracV
;icl3V=[14,5,8,6,30,7,21,7,17,13]
;icl3Verr=[5,3,6,7,14,2,9,3,6,5]
;icl3r=[14,15,10,8,31,7,17,12,16,16]
;icl3rerr=[7,6,3,7,12,3,9,5,3,5]

;ratios of elliptical profiles at the redshifts of clusters, ie no
;correction factors, and their errors
ratior = 100.*(lumr/ (lumr+clusterrflux))
ratiorerr = (ratior*rfracerr)


ratioB = 100.*(lumB / (lumB + clusterBflux))
ratioBerr = (ratioB*Bfracerr)

;print, "ratios",ratior, ratioB
print, "ratio errss",ratiorerr, ratioBerr

;color, colorerror, and color gradient
color=[1.7,2.4,2.3,2.1,2.15,0.5,-.1,0.95,0.2,0.8]; B-r or V-r, no corrections at all
color =[1.89,2.63,2.54,2.48,2.54,1.97,1.49,2.72,2.15,2.75]  ; k+ ev + v->B corrected to z=0
colorerr = [0.3,0.5,0.6,0.3,0.3,0.2,0.3,0.3,0.7,0.1]
;colorgrad = [1,2,2,0,2,2,0,-1,0,-1]
colorgrad = [1,1,1,1,2,2,0,-1,1,-1]
;non weighted error array
err = findgen(n_elements(colorgrad)) -  findgen(n_elements(colorgrad)) + 1.0

;-------------------------------------------------------------------
;plot, cd, color, psym = 2, thick = 3, xtitle = "cD? 1=yes, 0=no", $
;  ytitle = "B-r",xrange = [-0.5,1.5],yrange=[1.0,3.0],ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;
;;oplot, newx,newcolory, thick = 3;, psym = 6       
;errplot, cd, color-colorerr, color+colorerr
;spearman = r_correlate(cd, color)
;xyouts, 1.01, 2.8, string(spearman(0)), charthick = 3
;
;start=[-1E11,7.]
;result = MPFITFUN('linear',cd, color,colorerr, start,perror =test,/quiet)
;print, test
;oplot, cd,result(0)*cd + result(1), color = colors.black
;oplot, cd, (result(0) -2*test(0))*cd + (result(1) - 2*test(1)), color = colors.black
;oplot, cd, (result(0)+2*test(0) )*cd + (result(1) + 2*test(1)), color = colors.black
;;
;plot, cd, iclfracV, psym = 3, thick = 3, xtitle = "cD? 1=yes, 0=no", $
;  ytitle = "icl fraction (%) ",xrange = [-0.5,1.5],yrange = [0,35],$
;   charthick = 3, $
;  xthick = 3, ythick = 3;
;
;oplot, cd[5:9], iclfracV[5:9], color = colors.green, psym = 2, thick =3
;oplot, cd[0:4], iclfracV[0:4], color = colors.blue, psym = 2, thick =3
;errplot, cd,iclfracV-iclfracVerr, iclfracV+iclfracVerr
;
;oplot, cd, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, cd,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;
;oplot, newx,newfracVy, thick = 3;, psym = 6       
;oplot, newx,newfracry, thick = 3;, psym = 6       
;spearman = r_correlate(cd, iclfracr)
;xyouts, 1.01, 34, string(spearman(0)), charthick = 3
;
;start=[-1E11,7E11]
;result = MPFITFUN('linear',cd, iclfracr,iclfracrerr, start,perror =test,/quiet)
;print, test
;oplot, cd,result(0)*cd + result(1), color = colors.black
;oplot, cd, (result(0) -2*test(0))*cd + (result(1) - 2*test(1)), color = colors.black
;oplot, cd, (result(0)+2*test(0) )*cd + (result(1) + 2*test(1)), color = colors.black
;;;
plot, cd, ratioB, psym = 3, thick = 3, xtitle = "cD? 1=yes, 0=no", $
  ytitle = "icl 0.3 fraction (%) ",xrange = [-0.5,1.5],yrange = [0,30],$
   charthick = 3, $
  xthick = 3, ythick = 3
oplot, cd[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
oplot, cd[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
errplot, cd,ratioB-ratioBerr, ratioB+ratioBerr;

oplot, cd, ratior, psym = 4, thick = 3, color = colors.black
errplot, cd,ratior-ratiorerr, ratior+ratiorerr
spearman = r_correlate(cd, ratior)
xyouts, 1.01, 26, string(spearman(0), format='(F10.2)'), charthick = 3

;oplot, newx,newfrac3Vy, thick = 3;, psym = 6       
;oplot, newx,newfrac3ry, thick = 3;, psym = 6       

start=[-1E11,7E11]
result = MPFITFUN('linear',cd, ratior,ratiorerr, start,perror =test,/quiet)
print, test
oplot, cd,result(0)*cd + result(1), color = colors.black
oplot, cd, (result(0) -2*test(0))*cd + (result(1) - 2*test(1)), color = colors.black
oplot, cd, (result(0)+2*test(0) )*cd + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, cd, colorgrad, psym = 2, thick = 3, xtitle = "cD? 1=yes, 0=no", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],xrange=[-0.5,1.5],$
;  xthick = 3, ystyle = 1,ythick = 3;
;spearman = r_correlate(cd, colorgrad)
;xyouts, 1.01, 2.0, string(spearman(0)), charthick = 3
;start=[-1E11,7.]
;result = MPFITFUN('linear',cd, colorgrad,err, start,perror =test,/quiet)
;;print, "cd colorgrad",test
;oplot, cd,result(0)*cd + result(1), color = colors.black
;oplot, cd, (result(0) -2*test(0))*cd + (result(1) - 2*test(1)), color = colors.black
;oplot, cd, (result(0)+2*test(0) )*cd + (result(1) + 2*test(1)), color = colors.black;
;
;oplot, newx, newcolorgrady, thick = 3
;;
;plot, cd, lumr, psym = 2, thick = 3, xtitle = "cD? 1=yes, 0=no", $
;  ytitle = "ICL r-band luminosity",xrange = [-0.5,1.5],yrange = [1E10,23E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, cd, lumr-lumrerr, lumr+lumrerr
;spearman = r_correlate(cd, lumr)
;xyouts, 1.01, 1.3E12, string(spearman(0)), charthick = 3
;
;start=[-1E11,7E11]
;result = MPFITFUN('linear',cd, lumr,lumrerr, start,perror =test,/quiet)
;print, test
;oplot, cd,result(0)*cd + result(1), color = colors.black
;oplot, cd, (result(0) -2*test(0))*cd + (result(1) - 2*test(1)), color = colors.black
;oplot, cd, (result(0)+2*test(0) )*cd + (result(1) + 2*test(1)), color = colors.black

;;;
plot, cd, lumrzke, psym = 2, thick = 3, xtitle = "cD? 1=yes, 0=no", $
  ytitle = "ICL r-band luminosity, flux+k+ev corr",xrange = [-0.5,1.5],yrange = [0,10E11],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
errplot, cd, lumrzke-lumrzkeerr, lumrzke+lumrzkeerr
spearman = r_correlate(cd, lumrzke)
xyouts, 1.01, 9E11, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',cd, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;print, test
oplot, cd,result(0)*cd + result(1), color = colors.black
oplot, cd, (result(0) -2*test(0))*cd + (result(1) - 2*test(1)), color = colors.black
oplot, cd, (result(0)+2*test(0) )*cd + (result(1) + 2*test(1)), color = colors.black

;oplot, newx, newlumy, thick = 3
;--------------------------------------------------------------------

;plot, mass, color, psym = 2, thick = 3, xtitle = "mass h70 E14", $
;  ytitle = "B-r",yrange=[1.0,3.0], ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oploterror, mass, color , masserr, colorerr, psym = 3, errthick = 2
;spearman = r_correlate(mass, color)
;xyouts, 32, 2.8, string(spearman(0)), charthick = 3;;
;
;start=[-1E11,7.]
;result = MPFITFUN('linear',mass, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, mass,result(0)*mass + result(1), color = colors.black
;oplot, mass, (result(0) -2*test(0))*mass + (result(1) - 2*test(1)), color = colors.black
;oplot, mass, (result(0)+2*test(0) )*mass + (result(1) + 2*test(1)), color = colors.black
;;;;
;plot, mass, iclfracV, psym = 3, thick = 3, xtitle = "mass h70 E14", $
;  ytitle = "icl fraction (%) ",yrange = [0,70],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, mass[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, mass[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue;
;
;oploterror, mass, iclfracV, masserr, iclfracVerr, psym = 3, errthick = 2
;
;oplot, mass, iclfracr, psym = 4, thick = 3, color = colors.black
;oploterror, mass, iclfracr, masserr, iclfracrerr, psym = 3, errthick = 2
;
linx = [1,2,3,4,5,10,11,40]
linx = linx *1E14
liny =[30,40,45,49,52,64,65,75]
murantex = [1,2,3,4,5,7,40]
murantex = murantex*1E14
murantey=[25.5,33,37,39,41,45,50]
;oplot, linx,liny, thick = 3
;oplot, murantex,murantey, thick=3, linestyle=2
;xyouts, 35, 75,"lin04",charthick=3
;xyouts, 34, 51,"murante04",charthick=3
;spearman = r_correlate(mass, iclfracr)
;xyouts, 32, 34, string(spearman(0)), charthick = 3
;;;
plot, mass, ratioB, psym = 3, thick = 3, xtitle = "Cluster Mass (solar masses)", $
  ytitle = "ICL fraction (%) ", yrange = [0,80],$
   charthick = 3,$
  xthick = 3, ythick = 3
oplot, mass[5:9], ratioB[5:9], psym = 4, thick = 3;, color = colors.green
oplot, mass[0:4], ratioB[0:4], psym = 6, thick = 3;, color = colors.blue
oploterror, mass, ratioB, masserr, ratioBerr, psym = 3, errthick = 2

oplot, mass, ratior, psym = 2, thick = 3;, color = colors.black
oploterror, mass, ratior, masserr, ratiorerr, psym = 3, errthick = 2
spearman = r_correlate(mass, ratior)
;xyouts, 32, 26, string(spearman(0)), charthick = 3
oplot, linx,liny, thick = 3
oplot, murantex,murantey, thick=3, linestyle=2
xyouts, 30E14, 75,"Lin & Mohr 04",charthick=3
xyouts, 29E14, 51,"Murante et al. 04",charthick=3

xyouts, 11E14, 65, "X", alignment=0.5,charthick=3
xyouts, 7E14, 45, "X", alignment=0.5,charthick = 3

;;;
;plot, mass, colorgrad, psym = 2, thick = 3, xtitle = "mass h70 E14", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],$
;  xthick = 3, ystyle = 1,ythick = 3
;spearman = r_correlate(mass, colorgrad)
;xyouts, 32, 2.0, string(spearman(0)), charthick = 3;

;start=[-1E11,7.]
;result = MPFITFUN('linear',mass, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, mass,result(0)*mass + result(1), color = colors.black
;oplot, mass, (result(0) -2*test(0))*mass + (result(1) - 2*test(1)), color = colors.black
;oplot, mass, (result(0)+2*test(0) )*mass + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, mass, lumr, psym = 2,thick = 3, xtitle = "mass h70 E14", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,15E11],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oploterror, mass, lumr, masserr, lumrzkeerr, psym = 3, errthick=2
;spearman = r_correlate(mass, lumr)
;xyouts, 32, 9E11, string(spearman(0)), charthick = 3;

;start=[-1E11,7E11]
;result = MPFITFUN('linear',mass, lumr,lumrzkeerr, start,perror =test,/quiet)
;print, test
;oplot, mass,result(0)*mass + result(1), color = colors.black
;oplot, mass, (result(0) -2*test(0))*mass + (result(1) - 2*test(1)), color = colors.black
;oplot, mass, (result(0)+2*test(0) )*mass + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, mass, lumrzke, psym = 2,thick = 3, xtitle = "mass h70 E14", $
;  ytitle = "ICL r-band luminosity flux+k+ev corr",yrange = [0,10E11],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oploterror, mass, lumrzke, masserr, lumrzkeerr, psym = 3, errthick=2
;spearman = r_correlate(mass, lumrzke)
;xyouts, 28, 9E11, string(spearman(0)), charthick = 3
;
;start=[-1E11,7E11]
;result = MPFITFUN('linear',mass, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;;print, test
;oplot, mass,result(0)*mass + result(1), color = colors.black
;oplot, mass, (result(0) -2*test(0))*mass + (result(1) - 2*test(1)), color = colors.black
;oplot, mass, (result(0)+2*test(0) )*mass + (result(1) + 2*test(1)), color = colors.black

;;;;
;;--------------------------------------------------------------------

plot, z, color, psym = 2, thick = 3, xtitle = "Redshift", $
  ytitle = "B-r (magnitudes)",yrange = [1.0,3.0], ystyle = 1,$
   charthick = 3, $
  xthick = 3, ythick = 3
errplot, z, color + colorerr, color-colorerr
spearman = r_correlate(z, color)
;xyouts, 0.28, 2.8, string(spearman(0)), charthick = 3


filelist10 = [ '/Users/jkrick/umich/icl/color/m22z12.color',$
              '/Users/jkrick/umich/icl/color/m42z12.color',$
              '/Users/jkrick/umich/icl/color/m22z8.color',$
              '/Users/jkrick/umich/icl/color/m42z8.color',$
              '/Users/jkrick/umich/icl/color/m72z13.color']

FOR n = 0, 4  DO begin
    OPENR, ssplun,filelist10[n], /GET_LUN

;read in the radial profile
    rows= 50
    zarr= FLTARR(rows)
    agearr= FLTARR(rows)
    cevarr= FLTARR(rows)
    crfarr= FLTARR(rows)
    cnearr = FLTARR(rows)
    karr= FLTARR(rows)
    ekarr= FLTARR(rows)
    
    zed = 0.0			
    age = 0.0
    cev = 0.0
    junk = 0.0
    
    j = 0
    FOR j=0,rows - 1 DO BEGIN
        READF, ssplun, zed,junk,age,junk,crf,cne,cev,junk,junk,junk,ek,k
        zarr(j) = zed
        agearr(j) = age
        cevarr(j) = cev - 0.3
        cnearr(j) = cne
        crfarr(j) = crf
        karr(j) = k
        ekarr(j) = ek
    ENDFOR
    
    close, ssplun
    free_lun, ssplun
    
    
;    oplot,zarr,cevarr,color = colors.magenta, thick = 3, linestyle = n

ENDFOR



;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[-1E11,7.]
result = MPFITFUN('linear',z, color,colorerr, start,perror =test,/quiet)
print, "color vs redshift", test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, z, iclfracV, psym = 3, thick = 3, xtitle = "redshift", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, z[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, z[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, z,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, z, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, z,iclfracr-iclfracrerr, iclfracr+iclfracrerr
willmanx=[0.,0.05,0.1,0.2,0.3,0.35]
willmany =[19,19,17,16.5,13.9,13]


;oplot, willmanx, willmany, thick=3
;xyouts, 0.005,19.2,"Willman04", charthick=3
;spearman = r_correlate(z, iclfracr)
;xyouts, 0.28, 26, string(spearman(0)), charthick = 3
;;;
plot, z, ratioB, psym = 3, thick = 3, xtitle = "Redshift", $
  ytitle = "ICL fraction (%) ", yrange = [0,30],$
   charthick = 3, $
  xthick = 3, ythick = 3
oplot, z[5:9], ratioB[5:9], psym = 4, thick = 3;, color = colors.green
oplot, z[0:4], ratioB[0:4], psym = 6, thick = 3;, color = colors.blue
errplot, z,ratioB-ratioBerr, ratioB+ratioBerr

oplot, z, ratior, psym = 2, thick = 3;, color = colors.black
errplot, z,ratior-ratiorerr, ratior+ratiorerr
spearman = r_correlate(z, ratior)
xyouts, 0.28, 26, string(spearman(0), format='(F10.2)'), charthick = 3
oplot, willmanx, willmany, thick=3
xyouts, 0.005,19.2,"Willman04", charthick=3
xyouts, 0.005,18.0,"at r200", charthick=3

start=[-1E11,7E11]
result = MPFITFUN('linear',z, ratior,ratiorerr, start,perror =test,/quiet)
;;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, z, colorgrad, psym = 2, thick = 3, xtitle = "redshift", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],ystyle = 1,$
;  xthick = 3, ythick = 3;
;spearman = r_correlate(z, colorgrad)
;xyouts, 0.28, 2.0, string(spearman(0)), charthick = 3

;start=[-1E11,7.]
;result = MPFITFUN('linear',z, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, z, lumr, psym = 2, thick = 3, xtitle = "redshift", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,15E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, z, lumr-lumrzkeerr, lumr + lumrzkeerr
;spearman = r_correlate(z, lumr)
;xyouts, 0.28, 9E11, string(spearman(0)), charthick = 3

;start=[-1E11,7E11]
;result = MPFITFUN('linear',z, lumr,lumrzkeerr, start,perror =test,/quiet)
;;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black
;;;
plot, z, lumrzke, psym = 2, thick = 3, xtitle = "Redshift", $
  ytitle = "ICL r-band luminosity (solar luminosities)",yrange = [0,8E11],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
errplot, z, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(z, lumrzke)
;xyouts, 0.24, 9E11, string(spearman(0)), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',z, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black

;;;
plot, z, clusterrzke, psym = 2, thick = 3, xtitle = "Redshift", $
  ytitle = "Galaxy r-band luminosity (solar luminosities)",yrange = [0,0.6E13],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3

errplot, z, clusterrzke-clusterrzkeerr, clusterrzke + clusterrzkeerr
spearman = r_correlate(z, clusterrzke)
xyouts, 0.24, 5.1E12, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',z, clusterrzke,clusterrzkeerr, start,perror =test,/quiet)
;print, test
oplot, z,result(0)*z + result(1), color = colors.black
oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black



;;;
plot, z, m1m3, psym = 2, thick = 3, xtitle = "Redshift", $
  ytitle = "M3-M1 (magnitudes)", yrange=[0.1,1.2],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
;errplot, z, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(z, m1m3)
;xyouts, 0.24, 1.0, string(spearman(0)), charthick = 3

start=[-1E11,7]
result = MPFITFUN('linear',z, m1m3,err, start,perror =test,/quiet)
;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black

;;;

plot, z, mass, psym = 2, thick = 3, xtitle = "Redshift", $
  ytitle = "Cluster Mass (solar masses)",$
   charthick = 3,$
  xthick = 3, ythick = 3

errplot, z, mass - masserr, mass+masserr
spearman = r_correlate(z, mass)
;xyouts, 0.24, 5.1E14, string(spearman(0)), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',z, mass,masserr, start,perror =test,/quiet)
;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black


plot, z, density2, psym = 2, thick = 3, xtitle = "Redshift", $
  ytitle = "Number of Galaxies", yrange=[50,300],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
;errplot, z, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(z, m1m3)
;xyouts, 0.24, 1.0, string(spearman(0)), charthick = 3

start=[-1E11,7]
result = MPFITFUN('linear',z, m1m3,err, start,perror =test,/quiet)
;print, test
;oplot, z,result(0)*z + result(1), color = colors.black
;oplot, z, (result(0) -2*test(0))*z + (result(1) - 2*test(1)), color = colors.black
;oplot, z, (result(0)+2*test(0) )*z + (result(1) + 2*test(1)), color = colors.black

;;;

;

;-----------------------------------------------------------------------

;plot, rich, color, psym = 2, thick = 3, xtitle = "richness, 3 is richer than 2", $
;  ytitle = "B-r",xrange = [0.5,3.5],$
;   charthick = 3, yrange=[1.0,3.0],ystyle = 1,$
;  xthick = 3, ythick = 3
;;oplot, newrichx,newrichcolory, thick = 3;, psym = 6       
;spearman = r_correlate(rich, color)
;xyouts, 2.9, 2.8, string(spearman(0)), charthick = 3

;start=[-1E11,7.]
;result = MPFITFUN('linear',rich, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, rich,result(0)*rich + result(1), color = colors.black
;oplot, rich, (result(0) -2*test(0))*rich + (result(1) - 2*test(1)), color = colors.black
;oplot, rich, (result(0)+2*test(0) )*rich + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, rich, iclfracV, psym = 3, thick = 3, xtitle = "richness, 3 is richer than 2", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],xrange = [0.5,3.5],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, rich[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, rich[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, rich,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, rich, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, rich,iclfracr-iclfracrerr, iclfracr+iclfracrerr

;oplot, newrichx,newrichfracVy, thick = 3;, psym = 6       
;oplot, newrichx,newrichfracry, thick = 3;, psym = 6       
;spearman = r_correlate(rich, iclfracr)
;xyouts, 2.9, 26, string(spearman(0)), charthick = 3
;;;
;plot, rich, ratioB, psym = 3, thick = 3, xtitle = "richness, 3 is richer than 2", $
;  ytitle = "icl0.3 fraction (%) ", yrange = [0,30],xrange = [-1.0,3.5],$
;  charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, rich[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
;oplot, rich[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, rich,ratioB-ratioBerr, ratioB+ratioBerr;
;
;oplot, rich, ratior, psym = 4, thick = 3, color = colors.black
;errplot, rich,ratior-ratiorerr, ratior+ratiorerr
;
;;oplot, newrichx,newrichfrac3Vy, thick = 3;, psym = 6       
;;oplot, newrichx,newrichfrac3ry, thick = 3;, psym = 6       
;spearman = r_correlate(rich, ratior)
;xyouts, 2.9, 26, string(spearman(0)), charthick = 3

;start=[-1E11,7.]
;result = MPFITFUN('linear',rich, ratior,ratiorerr, start,perror =test,/quiet)
;;print, test
;oplot, rich,result(0)*rich + result(1), color = colors.black
;oplot, rich, (result(0) -2*test(0))*rich + (result(1) - 2*test(1)), color = colors.black
;oplot, rich, (result(0)+2*test(0) )*rich + (result(1) + 2*test(1)), color = colors.black
;;;;
;plot, rich, colorgrad, psym = 2, thick = 3, xtitle = "richness", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],xrange=[-0.5,3.5],$
;  xthick = 3, ystyle = 1,ythick = 3
;
;oplot, newrichx,newrichcolorgrady, thick = 3;, psym = 6       
;spearman = r_correlate(rich, colorgrad)
;xyouts, 2.9, 2.0, string(spearman(0)), charthick = 3

;start=[-1E11,7.]
;result = MPFITFUN('linear',rich, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, rich,result(0)*rich + result(1), color = colors.black
;oplot, rich, (result(0) -2*test(0))*rich + (result(1) - 2*test(1)), color = colors.black
;oplot, rich, (result(0)+2*test(0) )*rich + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, rich, lumr, psym = 2, thick = 3, xtitle = "richness", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,15E11],$
;   charthick = 3, xrange=[-0.5,3.5],$
;  xthick = 3, ythick = 3

;oplot, newrichx,newrichlumy, thick = 3;, psym = 6       
;errplot, rich, lumr-lumrerr, lumr+lumrerr
;spearman = r_correlate(rich, lumr)
;xyouts, 2.9, 9E11, string(spearman(0)), charthick = 3

;start=[-1E11,7E11]
;result = MPFITFUN('linear',rich, lumr,lumrerr, start,perror =test,/quiet)
;;print, test
;oplot, rich,result(0)*rich + result(1), color = colors.black
;oplot, rich, (result(0) -2*test(0))*rich + (result(1) - 2*test(1)), color = colors.black
;oplot, rich, (result(0)+2*test(0) )*rich + (result(1) + 2*test(1)), color = colors.black
;;;
plot, rich, lumrzke, psym = 2, thick = 3, xtitle = "Richness Class", $
  ytitle = "ICL r-band Luminosity (solar luminosities)",yrange = [0,8E11],$
   charthick = 3, xrange=[-0.5,3.5],$
  xthick = 3, ythick = 3

;oplot, newrichx,newrichlumy, thick = 3;, psym = 6       
errplot, rich, lumrzke-lumrzkeerr, lumrzke+lumrzkeerr
spearman = r_correlate(rich, lumrzke)
xyouts, -0.5, 7E11, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',rich, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;print, test
oplot, rich,result(0)*rich + result(1), color = colors.black
oplot, rich, (result(0) -2*test(0))*rich + (result(1) - 2*test(1)), color = colors.black
oplot, rich, (result(0)+2*test(0) )*rich + (result(1) + 2*test(1)), color = colors.black

;-------------------------------------------------------------

;plot, sigma, mass, psym = 2, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "mass",xrange = [600,2000],$
;   charthick = 3, $
;  xthick = 3, ythick = 3

;oploterror, sigma, mass, sigmaerr, masserr, psym=3, errthick=2
;spearman = r_correlate(sigma, mass)
;xyouts, 1600, 32., string(spearman(0)), charthick = 3

;start=[1500.,7E12]
;result = MPFITFUN('linear',sigma, mass,masserr, start,perror =test,/quiet)
;;print, test
;oplot, sigma,result(0)*sigma + result(1), color = colors.black
;oplot, sigma, (result(0) -2*test(0))*sigma + (result(1) - 2*test(1)), color = colors.black
;oplot, sigma, (result(0)+2*test(0) )*sigma + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, sigma, iclfracV, psym = 3, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],xrange = [600,2000],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, sigma, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, sigma,iclfracr-iclfracrerr, iclfracr+iclfracrerr

;oplot, sigma[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, sigma[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, sigma,iclfracV-iclfracVerr, iclfracV+iclfracVerr
;spearman = r_correlate(sigma, iclfracr)
;xyouts, 1600, 26, string(spearman(0)), charthick = 3

;start=[1500.,7.]
;result = MPFITFUN('linear',sigma, iclfracr,iclfracrerr, start,perror =test,/quiet)
;print, "sigma iclfracr",test
;oplot, sigma,result(0)*sigma + result(1), color = colors.black
;oplot, sigma, (result(0) -2*test(0))*sigma + (result(1) - 2*test(1)), color = colors.black
;oplot, sigma, (result(0)+2*test(0) )*sigma + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, sigma, ratioB, psym = 3, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "icl0.3 fraction (%) ", yrange = [0,30],xrange = [600,2000],$
;  charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, sigma, ratior, psym = 4, thick = 3, color = colors.black
;errplot, sigma,ratior-ratiorerr, ratior+ratiorerr
;
;oplot, sigma[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
;oplot, sigma[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, sigma,ratioB-ratioBerr, ratioB+ratioBerr
;spearman = r_correlate(sigma, ratior)
;xyouts, 1600, 26, string(spearman(0)), charthick = 3
;
;start=[1500.,7.]
;result = MPFITFUN('linear',sigma, ratior,ratiorerr, start,perror =test,/quiet)
;print, "sigms ratior",test
;oplot, sigma,result(0)*sigma + result(1), color = colors.black
;oplot, sigma, (result(0) -2*test(0))*sigma + (result(1) - 2*test(1)), color = colors.black
;oplot, sigma, (result(0)+2*test(0) )*sigma + (result(1) + 2*test(1)), color = colors.black

;;;
;plot, sigma, color, psym = 2, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "B-r",xrange = [600,2000],$
;   charthick = 3, yrange=[1.0,3.0],ystyle=1,$
;  xthick = 3, ythick = 3
;oploterror, sigma, color, sigmaerr, colorerr, psym = 3, errthick=3
;spearman = r_correlate(sigma, color)
;xyouts, 1600, 2.8, string(spearman(0)), charthick = 3
;

;;;
;plot, sigma, colorgrad, psym = 2, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",xrange = [600,2000],$
;   charthick = 3, yrange=[-0.5,2.5], ystyle = 1,$
;  xthick = 3, ythick = 3
;junk = fltarr(10)
;oploterror, sigma, colorgrad, sigmaerr, junk,psym=3, errthick=2
;spearman = r_correlate(sigma, colorgrad)
;xyouts, 1600, 2.0, string(spearman(0)), charthick = 3

;start=[1500,7]
;result = MPFITFUN('linear',sigma, colorgrad,err, start,perror =test,/quiet)
;print, "sigma colorgrad", test
;oplot, sigma,result(0)*sigma + result(1), color = colors.black
;oplot, sigma, (result(0) -2*test(0))*sigma + (result(1) - 2*test(1)), color = colors.black
;oplot, sigma, (result(0)+2*test(0) )*sigma + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, sigma, lumr, psym = 2, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,15E11],$
;   charthick = 3,xrange = [600,2000], $
;  xthick = 3, ythick = 3;

;oploterror, sigma, lumr, sigmaerr, lumrzkeerr, psym=3, errthick=2
;spearman = r_correlate(sigma, lumr)
;xyouts, 1600, 1E12, string(spearman(0)), charthick = 3
;;;
;plot, sigma, lumrzke, psym = 2, thick = 3, xtitle = "velocity dispersion", $
;  ytitle = "ICL r-band luminosity flux+k+ev corr",yrange = [0,10E11],$
;   charthick = 3,xrange = [600,2000], $
;  xthick = 3, ythick = 3
;
;oploterror, sigma, lumrzke, sigmaerr, lumrzkeerr, psym=3, errthick=2
;spearman = r_correlate(sigma, lumrzke)
;xyouts, 1600, 9E11, string(spearman(0)), charthick = 3

;----------------------------------------------------------------------------

;plot, m1m3, color, psym = 2, thick = 3, xtitle = "M3-M1", $
;  ytitle = "B-r",yrange = [1.0,3.0], ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;errplot, m1m3, color + colorerr, color-colorerr
;spearman = r_correlate(m1m3, color)
;xyouts, 0.9, 2.8, string(spearman(0)), charthick = 3
;
;start=[-1E11,7]
;result = MPFITFUN('linear',m1m3, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
;oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, m1m3, iclfracV, psym = 3, thick = 3, xtitle = "M3-M1", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, m1m3[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, m1m3[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, m1m3,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, m1m3, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, m1m3,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;spearman = r_correlate(m1m3, iclfracr)
;xyouts, 0.9, 26, string(spearman(0)), charthick = 3
;;;
plot, m1m3, ratioB, psym = 3, thick = 3, xtitle = "M3-M1", $
  ytitle = "icl 0.3 fraction (%) ", yrange = [0,30],$
   charthick = 3, $
  xthick = 3, ythick = 3
oplot, m1m3[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
oplot, m1m3[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
errplot, m1m3,ratioB-ratioBerr, ratioB+ratioBerr;

oplot, m1m3, ratior, psym = 4, thick = 3, color = colors.black
errplot, m1m3,ratior-ratiorerr, ratior+ratiorerr

spearman = r_correlate(m1m3, ratior)
xyouts, 0.9, 26, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7]
result = MPFITFUN('linear',m1m3, ratior,ratiorerr, start,perror =test,/quiet)
;;print, test
oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black
;;;;
;;plot, m1m3, colorgrad, psym = 2, thick = 3, xtitle = "M3-M1", $
;;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;;   charthick = 3, yrange=[-0.5,2.5],ystyle = 1,$
;;  xthick = 3, ythick = 3;
;spearman = r_correlate(m1m3, colorgrad)
;xyouts, 0.9, 2.0, string(spearman(0)), charthick = 3
;;;
;plot, m1m3, lumr, psym = 2, thick = 3, xtitle = "M3-M1", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,15E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, m1m3, lumr-lumrerr, lumr + lumrerr
;spearman = r_correlate(m1m3, lumr)
;xyouts, 0.9, 1E12, string(spearman(0)), charthick = 3

;start=[-1E11,7E11]
;result = MPFITFUN('linear',m1m3, lumr,lumrerr, start,perror =test,/quiet)
;;print, test
;oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
;oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black
;;;
;coeff1 = ROBUST_LINEFIT(m1m3,lumr, yfit, sig, coeff_sig)
;oplot, m1m3,yfit, thick =3
;print, "sigma", sig
;oplot, m1m3,yfit-sig, linestyle=2
;oplot, m1m3,yfit+sig, linestyle=2
;start=[-1E11,7E11]
;result = MPFITFUN('linear',m1m3, lumr,lumrerr, start)
;oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black

plot, m1m3, lumrzke, psym = 2, thick = 3, xtitle = "M3-M1 (magnitudes)", $
  ytitle = "ICL r-band Luminosity (solar luminosities)",yrange = [0,8E11],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
errplot, m1m3, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(m1m3, lumrzke)
xyouts, 0.9, 7E11, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',m1m3, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;print, test
;oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
;oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black

;;;;;;;;

plot, m1m3, density2, psym = 2, thick = 3, xtitle = "M3-M1 (magnitudes)", $
  ytitle = "Number of Galaxies",yrange = [50,300],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
;errplot, m1m3, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(m1m3, density2)
xyouts, 0.9, 260, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1,50]
result = MPFITFUN('linear',m1m3, density2,err, start,perror =test)
;print, test
;oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
;oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black

;;;;;;;;

plot, m1m3, rich, psym = 2, thick = 3, xtitle = "M3-M1", $
  ytitle = "richness",yrange = [-1,4],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
;errplot, m1m3, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(m1m3, rich)
xyouts, 0.9, 200, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,1]
result = MPFITFUN('linear',m1m3, rich,err, start,perror =test,/quiet)
;print, test
oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black

;;;;;;;;

plot, m1m3, rvir, psym = 2, thick = 3, xtitle = "M3-M1", $
  ytitle = "one third rvir",$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
;errplot, m1m3, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(m1m3, rvir)
xyouts, 0.9, 200, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,1]
result = MPFITFUN('linear',m1m3, rvir,err, start,perror =test,/quiet)
;print, test
oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black

;;;;;;;;

plot, m1m3, clusterrzke, psym = 2, thick = 3, xtitle = "M3-M1 (magnitudes)", $
  ytitle = "Galaxy r-band luminosity (solar luminosities)",$
   charthick = 3, ystyle = 1, yrange=[0,5E12],$
  xthick = 3, ythick = 3
;errplot, m1m3, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(m1m3, clusterrzke)
xyouts, 0.9, 4E12, string(spearman(0), format='(F10.2)'), charthick = 3
errplot, m1m3, clusterrzke-clusterrzkeerr, clusterrzke + clusterrzkeerr

start=[-1E11,1E11]
result = MPFITFUN('linear',m1m3, clusterrzke ,clusterrzkeerr, start,perror =test,/quiet)
;print, test
;oplot, m1m3,result(0)*m1m3 + result(1), color = colors.black
;oplot, m1m3, (result(0) -2*test(0))*m1m3 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m3, (result(0)+2*test(0) )*m1m3 + (result(1) + 2*test(1)), color = colors.black


;----------------------------------------------------------------------------

;plot, m1m2, color, psym = 2, thick = 3, xtitle = "M2-M1", $
;  ytitle = "B-r",yrange = [1.0,3.0], ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;errplot, m1m2, color + colorerr, color-colorerr
;spearman = r_correlate(m1m2, color)
;xyouts, 0.75, 2.8, string(spearman(0)), charthick = 3

;start=[-1E11,7]
;result = MPFITFUN('linear',m1m2, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, m1m2,result(0)*m1m2 + result(1), color = colors.black
;oplot, m1m2, (result(0) -2*test(0))*m1m2 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m2, (result(0)+2*test(0) )*m1m2 + (result(1) + 2*test(1)), color = colors.black

;;;
;plot, m1m2, iclfracV, psym = 3, thick = 3, xtitle = "M2-M1", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, m1m2[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, m1m2[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, m1m2,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, m1m2, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, m1m2,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;spearman = r_correlate(m1m2, iclfracr)
;xyouts, 0.75, 26, string(spearman(0)), charthick = 3
;;;
;plot, m1m2, ratioB, psym = 3, thick = 3, xtitle = "M2-M1", $
;  ytitle = "icl 0.3 fraction (%) ", yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, m1m2[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
;oplot, m1m2[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, m1m2,ratioB-ratioBerr, ratioB+ratioBerr

;oplot, m1m2, ratior, psym = 4, thick = 3, color = colors.black
;errplot, m1m2,ratior-ratiorerr, ratior+ratiorerr
;spearman = r_correlate(m1m2, ratior)
;xyouts, 0.75, 26, string(spearman(0)), charthick = 3

;start=[-0.1,7]
;result = MPFITFUN('linear',m1m2, ratior,ratiorerr, start,perror =test,/quiet)
;;print, test
;oplot, m1m2,result(0)*m1m2 + result(1), color = colors.black;
;oplot, m1m2, (result(0) -2*test(0))*m1m2 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m2, (result(0)+2*test(0) )*m1m2 + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, m1m2, colorgrad, psym = 2, thick = 3, xtitle = "M2-M1", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],ystyle = 1,$
;  xthick = 3, ythick = 3;
;spearman = r_correlate(m1m2, colorgrad)
;xyouts, 0.75, 2.0, string(spearman(0)), charthick = 3
;;;
;plot, m1m2, lumr, psym = 2, thick = 3, xtitle = "M2-M1", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,15E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, m1m2, lumr-lumrerr, lumr + lumrerr
;spearman = r_correlate(m1m2, lumr)
;xyouts, 0.75, 1E12, string(spearman(0)), charthick = 3

;start=[-1E11,7E11]
;result = MPFITFUN('linear',m1m2, lumr,lumrerr, start,perror =test,/quiet)
;;print, test
;oplot, m1m2,result(0)*m1m2 + result(1), color = colors.black
;oplot, m1m2, (result(0) -2*test(0))*m1m2 + (result(1) - 2*test(1)), color = colors.black
;oplot, m1m2, (result(0)+2*test(0) )*m1m2 + (result(1) + 2*test(1)), color = colors.black

;;;
;coeff1 = ROBUST_LINEFIT(m1m2,lumr, yfit, sig, coeff_sig)
;oplot, m1m2,yfit, thick =3
;print, "sigma", sig
;oplot, m1m2,yfit-sig, linestyle=2
;oplot, m1m2,yfit+sig, linestyle=2
;start=[-1E11,7E11]
;result = MPFITFUN('linear',m1m2, lumr,lumrerr, start)
;oplot, m1m2,result(0)*m1m2 + result(1), color = colors.black

plot, m1m2, lumrzke, psym = 2, thick = 3, xtitle = "M2-M1", $
  ytitle = "ICL r-band luminosity flux+k+ev corr",yrange = [0,8E11],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
errplot, m1m2, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
spearman = r_correlate(m1m2, lumrzke)
xyouts, 0.75, 9E11, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',m1m2, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;print, test
oplot, m1m2,result(0)*m1m2 + result(1), color = colors.black
oplot, m1m2, (result(0) -2*test(0))*m1m2 + (result(1) - 2*test(1)), color = colors.black
oplot, m1m2, (result(0)+2*test(0) )*m1m2 + (result(1) + 2*test(1)), color = colors.black;


;----------------------------------------------------------------------------


;plot, clusterrzke, color, psym = 2, thick = 3, xtitle = "galaxy r flux", $
;  ytitle = "B-r",yrange = [1.0,3.0], ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;;errplot, clusterrzke, color + colorerr, color-colorerr
;oploterror, clusterrzke, color , clusterrzkeerr, colorerr, psym = 3, errthick = 2
;spearman = r_correlate(clusterrzke, color)
;xyouts, 4.2E12, 2.8, string(spearman(0)), charthick = 3
;start=[-1E11,7]
;result = MPFITFUN('linear',clusterrzke, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, clusterrzke,result(0)*clusterrzke + result(1), color = colors.black
;oplot, clusterrzke, (result(0) -2*test(0))*clusterrzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterrzke, (result(0)+2*test(0) )*clusterrzke + (result(1) + 2*test(1)), color = colors.black
;;
;plot, clusterrzke, iclfracV, psym = 3, thick = 3, xtitle = "galaxy r flux", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, clusterrzke[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, clusterrzke[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, clusterrzke,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, clusterrzke, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, clusterrzke,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;spearman = r_correlate(clusterrzke, iclfracr)
;xyouts, 4.9E11, 26, string(spearman(0)), charthick = 3
;;
plot, clusterrzke, ratioB, psym = 3, thick = 3, xtitle = "Galaxy r-band luminosity (solar luminosities)", $
  ytitle = "ICL Fraction (%) ", yrange = [0,30],$
   charthick = 3, $
  xthick = 3, ythick = 3
oplot, clusterrzke[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
oplot, clusterrzke[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, clusterrzke,ratioB-ratioBerr, ratioB+ratioBerr
oploterror, clusterrzke, ratioB , clusterrzkeerr, ratioBerr, psym = 3, errthick = 2

oplot, clusterrzke, ratior, psym = 4, thick = 3, color = colors.black
;errplot, clusterrzke,ratior-ratiorerr, ratior+ratiorerr
oploterror, clusterrzke, ratior , clusterrzkeerr, ratiorerr, psym = 3, errthick = 2
spearman = r_correlate(clusterrzke, ratior)
;xyouts, 4.2E12, 26, string(spearman(0)), charthick = 3

start=[-1E11,7]
result = MPFITFUN('linear',clusterrzke, ratior,ratiorerr, start,perror =test,/quiet)
;print, "colorgrad",test
;oplot, clusterrzke,result(0)*clusterrzke + result(1), color = colors.black
;oplot, clusterrzke, (result(0) -2*test(0))*clusterrzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterrzke, (result(0)+2*test(0) )*clusterrzke + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, clusterrzke, colorgrad, psym = 2, thick = 3, xtitle = "galaxy r flux", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],ystyle = 1,$
;  xthick = 3, ythick = 3;
;spearman = r_correlate(clusterrzke, colorgrad)
;xyouts, 4.2E12, 2.2, string(spearman(0)), charthick = 3;

;start=[-1E11,7]
;result = MPFITFUN('linear',clusterrzke, colorgrad,err, start,perror =test,/quiet)
;print, "colorgrad",test
;oplot, clusterrzke,result(0)*clusterrzke + result(1), color = colors.black
;oplot, clusterrzke, (result(0) -2*test(0))*clusterrzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterrzke, (result(0)+2*test(0) )*clusterrzke + (result(1) + 2*test(1)), color = colors.black

;;;
;plot, clusterrzke, lumr, psym = 2, thick = 3, xtitle = "galaxy r flux", $
;  ytitle = "ICL r-band luminosity",yrange = [-1E11,15E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, clusterrzke, lumr-lumrerr, lumr + lumrerr
;coeff1 = ROBUST_LINEFIT(clusterrzke,lumr, yfit, sig, coeff_sig)
;oplot, clusterrzke,yfit, thick =3
;print, "sigma", sig
;oplot, clusterrzke,yfit-sig, linestyle=2
;oplot, clusterrzke,yfit+sig, linestyle=2

;start=[-1E11,7E11]
;result = MPFITFUN('linear',clusterrzke, lumr,lumrerr, start,perror =test,/quiet)
;;print, test
;oplot, clusterrzke,result(0)*clusterrzke + result(1), color = colors.black
;oplot, clusterrzke, (result(0) -2*test(0))*clusterrzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterrzke, (result(0)+2*test(0) )*clusterrzke + (result(1) + 2*test(1)), color = colors.black

;spearman = r_correlate(clusterrzke, lumr)
;xyouts, 4.1.3E12, 1E12, string(spearman(0)), charthick = 3
;;;;
plot, clusterrzke, lumrzke, psym = 2, thick = 3, xtitle = "Galaxy r-band luminosity (solar luminosities)", $
  ytitle = "ICL r-band luminosity (solar luminosities)",yrange = [0,10E11],$
   charthick = 3, ystyle = 1,$
  xthick = 3, ythick = 3
;errplot, clusterrzke, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
oploterror, clusterrzke, lumrzke , clusterrzkeerr, lumrzkeerr, psym = 3, errthick = 2

spearman = r_correlate(clusterrzke, lumrzke)
xyouts, 3E12, 9E11, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
print, "cluster flux vs ICL flux"
result = MPFITFUN('linear',clusterrzke, lumrzke,lumrzkeerr, start,perror =test)
;print, test
oplot, clusterrzke,result(0)*clusterrzke + result(1), color = colors.black
oplot, clusterrzke, (result(0) -2*test(0))*clusterrzke + (result(1) - 2*test(1)), color = colors.black
oplot, clusterrzke, (result(0)+2*test(0) )*clusterrzke + (result(1) + 2*test(1)), color = colors.black

;----------------------------------------------------------------------------

;plot, clusterBzke, color, psym = 2, thick = 3, xtitle = "galaxy B flux", $
;  ytitle = "B-r",yrange = [1.0,3.0], ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;errplot, clusterBzke, color + colorerr, color-colorerr
;oploterror, clusterBzke, color , clusterBzkeerr, colorerr, psym = 3, errthick = 2

;spearman = r_correlate(clusterBzke, color)
;xyouts, 5.7E11, 2.8, string(spearman(0)), charthick = 3

;start=[-1E11,2.5]
;result = MPFITFUN('linear',clusterBzke, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, clusterBzke,result(0)*clusterBzke + result(1), color = colors.black
;oplot, clusterBzke, (result(0) -2*test(0))*clusterBzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterBzke, (result(0)+2*test(0) )*clusterBzke + (result(1) + 2*test(1)), color = colors.black

;;;
;plot, clusterBzke, iclfracV, psym = 3, thick = 3, xtitle = "galaxy B flux", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, clusterBzke[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, clusterBzke[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, clusterBzke,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, clusterBzke, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, clusterBzke,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;spearman = r_correlate(clusterBzke, iclfracr)
;xyouts, 5.7E11, 26, string(spearman(0)), charthick = 3

;;;;
;plot, clusterBzke, ratioB, psym = 3, thick = 3, xtitle = "galaxy B flux", $
;  ytitle = "icl 0.3 fraction (%) ", yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, clusterBzke[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
;oplot, clusterBzke[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;;;errplot, clusterBzke,ratioB-ratioBerr, ratioB+ratioBerr
;;oploterror, clusterBzke, ratioB , clusterBzkeerr, ratioBerr, psym = 3, errthick = 2

;oplot, clusterBzke, ratior, psym = 4, thick = 3, color = colors.black
;;errplot, clusterBzke,ratior-ratiorerr, ratior+ratiorerr
;oploterror, clusterBzke, ratior , clusterBzkeerr, ratiorerr, psym = 3, errthick = 2

;spearman = r_correlate(clusterBzke, ratior)
;xyouts, 5.7E11, 26, string(spearman(0)), charthick = 3

;start=[-1E11,20]
;result = MPFITFUN('linear',clusterBzke, ratior,ratiorerr, start,perror =test,/quiet)
;;print, test
;oplot, clusterBzke,result(0)*clusterBzke + result(1), color = colors.black
;oplot, clusterBzke, (result(0) -2*test(0))*clusterBzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterBzke, (result(0)+2*test(0) )*clusterBzke + (result(1) + 2*test(1)), color = colors.black

;;;;
;plot, clusterBzke, colorgrad, psym = 2, thick = 3, xtitle = "galaxy B flux", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5],ystyle = 1,$
;  xthick = 3, ythick = 3;
;;oploterror, clusterBzke, colorgrad , clusterBzkeerr, colorerr, psym = 3, errthick = 2

;spearman = r_correlate(clusterBzke, colorgrad)
;xyouts, 5.7E11, 2.2, string(spearman(0)), charthick = 3

;start=[-1E11,1.5]
;result = MPFITFUN('linear',clusterBzke, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, clusterBzke,result(0)*clusterBzke + result(1), color = colors.black
;oplot, clusterBzke, (result(0) -2*test(0))*clusterBzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterBzke, (result(0)+2*test(0) )*clusterBzke + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, clusterBzke, lumr, psym = 2, thick = 3, xtitle = "galaxy B flux", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,11E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, clusterBzke, lumr-lumrzkeerr, lumr + lumrerr

;coeff1 = ROBUST_LINEFIT(clusterBzke,lumr, yfit, sig, coeff_sig)
;oplot, clusterBzke,yfit, thick =3
;print, "sigma", sig
;oplot, clusterBzke,yfit-sig, linestyle=2
;oplot, clusterBzke,yfit+sig, linestyle=2

;start=[-1E11,7E11]
;result = MPFITFUN('linear',clusterBzke, lumr,lumrerr, start,perror =test,/quiet)
;;print, test
;oplot, clusterBzke,result(0)*clusterBzke + result(1), color = colors.black
;oplot, clusterBzke, (result(0) -2*test(0))*clusterBzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterBzke, (result(0)+2*test(0) )*clusterBzke + (result(1) + 2*test(1)), color = colors.black

;start=[4E11]
;result = MPFITFUN('noslope',clusterBzke, lumr,lumrerr, start,perror =test)
;;print, test
;oplot, clusterBzke,findgen(n_elements(clusterBzke)) - findgen(n_elements(clusterBzke)) + result(0), color = colors.blue
;oplot, clusterBzke, (result(0) -2*test(0))*clusterBzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterBzke, (result(0)+2*test(0) )*clusterBzke + (result(1) + 2*test(1)), color = colors.black

;spearman = r_correlate(clusterBzke, lumr)
;xyouts, 5.7E11, 1E12, string(spearman(0)), charthick = 3
;;;
;plot, clusterBzke, lumBzke, psym = 2, thick = 3, xtitle = "galaxy B flux", $
;  ytitle = "ICL B-band luminosity flux+k+ev corr",yrange = [0,1E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;;errplot, clusterBzke, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
;oploterror, clusterBzke, lumBzke , clusterBzkeerr, lumBzkeerr, psym = 3, errthick = 2;

;spearman = r_correlate(clusterBzke, lumBzke)
;xyouts, 5.7E11, .9E11, string(spearman(0)), charthick = 3

;start=[-1E11,7E11]
;result = MPFITFUN('linear',clusterBzke, lumBzke,lumBzkeerr, start,perror =test,/quiet)
;print, test
;oplot, clusterBzke,result(0)*clusterBzke + result(1), color = colors.black
;oplot, clusterBzke, (result(0) -2*test(0))*clusterBzke + (result(1) - 2*test(1)), color = colors.black
;oplot, clusterBzke, (result(0)+2*test(0) )*clusterBzke + (result(1) + 2*test(1)), color = colors.black

;-------------------------------------------------------


;plot, ngals, color, psym = 2, thick = 3, xtitle = "ngals inside 0.3rvir", $
;  ytitle = "B-r",yrange = [1.0,3.0], xrange=[50,250],xstyle=1,ystyle = 1,$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;errplot, ngals, color + colorerr, color-colorerr;

;spearman = r_correlate(ngals, color)
;xyouts, 150, 2.8, string(spearman(0)), charthick = 3

;start=[-1E11,2.5]
;result = MPFITFUN('linear',ngals, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, ngals,result(0)*ngals + result(1), color = colors.black
;oplot, ngals, (result(0) -2*test(0))*ngals + (result(1) - 2*test(1)), color = colors.black
;oplot, ngals, (result(0)+2*test(0) )*ngals + (result(1) + 2*test(1)), color = colors.black

;;;
;plot, ngals, iclfracV, psym = 3, thick = 3, xtitle = "ngals inside 0.3rvir", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, ngals[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, ngals[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, ngals,iclfracV-iclfracVerr, iclfracV+iclfracVerr

;oplot, ngals, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, ngals,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;spearman = r_correlate(ngals, iclfracr)
;xyouts, 150, 26, string(spearman(0)), charthick = 3

;;;;
;plot, ngals, ratioB, psym = 3, thick = 3, xtitle = "ngals inside 0.3rvir", $
;  ytitle = "icl 0.3 fraction (%) ", yrange = [0,30],$
;   charthick = 3, xrange=[50,250],xstyle=1, $
;  xthick = 3, ythick = 3
;oplot, ngals[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
;oplot, ngals[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, ngals,ratioB-ratioBerr, ratioB+ratioBerr
;
;oplot, ngals, ratior, psym = 4, thick = 3, color = colors.black
;errplot, ngals,ratior-ratiorerr, ratior+ratiorerr
;spearman = r_correlate(ngals, ratior)
;xyouts, 150, 26, string(spearman(0)), charthick = 3
;
;start=[-1E11,20]
;result = MPFITFUN('linear',ngals, ratior,ratiorerr, start,perror =test,/quiet)
;;;print, test
;oplot, ngals,result(0)*ngals + result(1), color = colors.black
;oplot, ngals, (result(0) -2*test(0))*ngals + (result(1) - 2*test(1)), color = colors.black
;oplot, ngals, (result(0)+2*test(0) )*ngals + (result(1) + 2*test(1)), color = colors.black
;
;;;;;
;plot, ngals, colorgrad, psym = 2, thick = 3, xtitle = "ngals inside 0.3rvir", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",$
;   charthick = 3, yrange=[-0.5,2.5], xrange=[50,250],xstyle=1,ystyle = 1,$
;  xthick = 3, ythick = 3;
;spearman = r_correlate(ngals, colorgrad)
;xyouts, 150, 2.2, string(spearman(0)), charthick = 3
;
;start=[-1E11,1.5]
;result = MPFITFUN('linear',ngals, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, ngals,result(0)*ngals + result(1), color = colors.black
;oplot, ngals, (result(0) -2*test(0))*ngals + (result(1) - 2*test(1)), color = colors.black
;oplot, ngals, (result(0)+2*test(0) )*ngals + (result(1) + 2*test(1)), color = colors.black
;;;;
;;plot, ngals, lumr, psym = 2, thick = 3, xtitle = "ngals inside 0.3rvir", $
;;  ytitle = "ICL r-band luminosity",yrange = [1E11,11E11],$
;;   charthick = 3, ystyle = 1,$
;;  xthick = 3, ythick = 3
;;errplot, ngals, lumr-lumrerr, lumr + lumrerr
;
;;coeff1 = ROBUST_LINEFIT(ngals,lumr, yfit, sig, coeff_sig)
;;oplot, ngals,yfit, thick =3
;;print, "sigma", sig
;;oplot, ngals,yfit-sig, linestyle=2
;;oplot, ngals,yfit+sig, linestyle=2
;
;;start=[-1E11,7E11]
;;result = MPFITFUN('linear',ngals, lumr,lumrerr, start,perror =test,/quiet)
;;;print, test
;;oplot, ngals,result(0)*ngals + result(1), color = colors.black
;;oplot, ngals, (result(0) -2*test(0))*ngals + (result(1) - 2*test(1)), color = colors.black
;;oplot, ngals, (result(0)+2*test(0) )*ngals + (result(1) + 2*test(1)), color = colors.black
;
;;start=[4E11]
;;result = MPFITFUN('noslope',ngals, lumr,lumrerr, start,perror =test)
;;;print, test
;;oplot, ngals,findgen(n_elements(ngals)) - findgen(n_elements(ngals)) + result(0), color = colors.blue
;;oplot, ngals, (result(0) -2*test(0))*ngals + (result(1) - 2*test(1)), color = colors.black
;;oplot, ngals, (result(0)+2*test(0) )*ngals + (result(1) + 2*test(1)), color = colors.black
;
;spearman = r_correlate(ngals, lumr)
;xyouts, 150, 1E12, string(spearman(0)), charthick = 3
;;;;
;plot, ngals, lumrzke, psym = 2, thick = 3, xtitle = "ngals inside 0.3rvir", $
;  ytitle = "ICL r-band luminosity flux+k+ev corr",yrange = [0,10E11],$
;   charthick = 3, xrange=[50,250],xstyle=1, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, ngals, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
;
;spearman = r_correlate(ngals, lumrzke)
;xyouts, 150, 9E11, string(spearman(0)), charthick = 3
;
;start=[-1E11,7E11]
;result = MPFITFUN('linear',ngals, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;;print, test
;oplot, ngals,result(0)*ngals + result(1), color = colors.black
;oplot, ngals, (result(0) -2*test(0))*ngals + (result(1) - 2*test(1)), color = colors.black
;oplot, ngals, (result(0)+2*test(0) )*ngals + (result(1) + 2*test(1)), color = colors.black
;
;;-------------------------------------------------------
;
;
;plot, density, color, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.3rvir", $
;  ytitle = "B-r",yrange = [1.0,3.0],xstyle=1,ystyle = 1,xrange=[20,65],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;errplot, density, color + colorerr, color-colorerr
;
;spearman = r_correlate(density, color)
;xyouts, 50, 2.8, string(spearman(0)), charthick = 3
;
;;start=[-1E11,2.5]
;;result = MPFITFUN('linear',density, color,colorerr, start,perror =test,/quiet)
;;;print, test
;;oplot, density,result(0)*density + result(1), color = colors.black
;;oplot, density, (result(0) -2*test(0))*density + (result(1) - 2*test(1)), color = colors.black
;;oplot, density, (result(0)+2*test(0) )*density + (result(1) + 2*test(1)), color = colors.black
;
;;;;
;;plot, density, iclfracV, psym = 3, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.3rvir", $
;;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;;   charthick = 3, $
;;  xthick = 3, ythick = 3
;;oplot, density[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;;oplot, density[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;;errplot, density,iclfracV-iclfracVerr, iclfracV+iclfracVerr
;
;;oplot, density, iclfracr, psym = 4, thick = 3, color = colors.black
;;errplot, density,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;;spearman = r_correlate(density, iclfracr)
;;xyouts, 150, 26, string(spearman(0)), charthick = 3
;
;;;;;
;plot, density, ratioB, psym = 3, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.3rvir", $
;  ytitle = "icl 0.3 fraction (%) ", yrange = [0,30],xrange=[20,65],$
;   charthick = 3,xstyle=1, $
;  xthick = 3, ythick = 3
;oplot, density[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
;oplot, density[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, density,ratioB-ratioBerr, ratioB+ratioBerr
;
;oplot, density, ratior, psym = 4, thick = 3, color = colors.black
;errplot, density,ratior-ratiorerr, ratior+ratiorerr
;spearman = r_correlate(density, ratior)
;xyouts, 50, 26, string(spearman(0)), charthick = 3
;
;start=[-1E11,20]
;result = MPFITFUN('linear',density, ratior,ratiorerr, start,perror =test,/quiet)
;;;print, test
;oplot, density,result(0)*density + result(1), color = colors.black
;oplot, density, (result(0) -2*test(0))*density + (result(1) - 2*test(1)), color = colors.black
;oplot, density, (result(0)+2*test(0) )*density + (result(1) + 2*test(1)), color = colors.black
;
;;;;;
;plot, density, colorgrad, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.3rvir", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",xrange=[20,65],$
;   charthick = 3, yrange=[-0.5,2.5],xstyle=1,ystyle = 1,$
;  xthick = 3, ythick = 3;
;spearman = r_correlate(density, colorgrad)
;xyouts, 50, 2.2, string(spearman(0)), charthick = 3
;
;start=[-1E11,1.5]
;result = MPFITFUN('linear',density, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, density,result(0)*density + result(1), color = colors.black
;oplot, density, (result(0) -2*test(0))*density + (result(1) - 2*test(1)), color = colors.black
;oplot, density, (result(0)+2*test(0) )*density + (result(1) + 2*test(1)), color = colors.black
;;;;
;;plot, density, lumr, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.3rvir", $
;;  ytitle = "ICL r-band luminosity",yrange = [1E11,11E11],$
;;   charthick = 3, ystyle = 1,$
;;  xthick = 3, ythick = 3
;;errplot, density, lumr-lumrerr, lumr + lumrerr
;
;;coeff1 = ROBUST_LINEFIT(density,lumr, yfit, sig, coeff_sig)
;;oplot, density,yfit, thick =3
;;print, "sigma", sig
;;oplot, density,yfit-sig, linestyle=2
;;oplot, density,yfit+sig, linestyle=2
;
;;start=[-1E11,7E11]
;;result = MPFITFUN('linear',density, lumr,lumrerr, start,perror =test,/quiet)
;;;print, test
;;oplot, density,result(0)*density + result(1), color = colors.black
;;oplot, density, (result(0) -2*test(0))*density + (result(1) - 2*test(1)), color = colors.black
;;oplot, density, (result(0)+2*test(0) )*density + (result(1) + 2*test(1)), color = colors.black
;
;;start=[4E11]
;;result = MPFITFUN('noslope',density, lumr,lumrerr, start,perror =test)
;;;print, test
;;oplot, density,findgen(n_elements(density)) - findgen(n_elements(density)) + result(0), color = colors.blue
;;oplot, density, (result(0) -2*test(0))*density + (result(1) - 2*test(1)), color = colors.black
;;oplot, density, (result(0)+2*test(0) )*density + (result(1) + 2*test(1)), color = colors.black
;
;spearman = r_correlate(density, lumr)
;xyouts, 50, 1E12, string(spearman(0)), charthick = 3
;;;;
;plot, density, lumrzke, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.3rvir", $
;  ytitle = "ICL r-band luminosity flux+k+ev corr",yrange = [0,10E11],$
;   charthick = 3,xstyle=1, ystyle = 1,xrange=[20,65],$
;  xthick = 3, ythick = 3
;errplot, density, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr
;
;spearman = r_correlate(density, lumrzke)
;xyouts, 50, 9E11, string(spearman(0)), charthick = 3
;
;start=[-1E11,7E11]
;result = MPFITFUN('linear',density, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;;print, test
;oplot, density,result(0)*density + result(1), color = colors.black
;oplot, density, (result(0) -2*test(0))*density + (result(1) - 2*test(1)), color = colors.black
;oplot, density, (result(0)+2*test(0) )*density + (result(1) + 2*test(1)), color = colors.black



;-------------------------------------------------------


;plot, density2, color, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.8Mpc", $
;  ytitle = "B-r",yrange = [1.0,3.0],xstyle=1,ystyle = 1,xrange=[50,160],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;errplot, density2, color + colorerr, color-colorerr
;
;spearman = r_correlate(density2, color)
;xyouts, 50, 2.8, string(spearman(0)), charthick = 3
;
;start=[-1E11,2.5]
;result = MPFITFUN('linear',density2, color,colorerr, start,perror =test,/quiet)
;;print, test
;oplot, density2,result(0)*density2 + result(1), color = colors.black
;oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
;oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black

;;
;plot, density2, iclfracV, psym = 3, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.8Mpc", $
;  ytitle = "icl fraction (%) ",yrange = [0,30],$
;   charthick = 3, $
;  xthick = 3, ythick = 3
;oplot, density2[5:9], iclfracV[5:9], psym = 2, thick = 3, color = colors.green
;oplot, density2[0:4], iclfracV[0:4], psym = 2, thick = 3, color = colors.blue
;errplot, density2,iclfracV-iclfracVerr, iclfracV+iclfracVerr
;
;oplot, density2, iclfracr, psym = 4, thick = 3, color = colors.black
;errplot, density2,iclfracr-iclfracrerr, iclfracr+iclfracrerr
;spearman = r_correlate(density2, iclfracr)
;xyouts, 150, 26, string(spearman(0)), charthick = 3

;;;
plot, density2, ratioB, psym = 3, thick = 3, xtitle = "Number of Galaxies", $
  ytitle = "ICL Fraction (%) ", yrange = [0,30],xrange=[50,300],$
   charthick = 3,xstyle=1, $
  xthick = 3, ythick = 3
oplot, density2[5:9], ratioB[5:9], psym = 2, thick = 3, color = colors.green
oplot, density2[0:4], ratioB[0:4], psym = 2, thick = 3, color = colors.blue
errplot, density2,ratioB-ratioBerr, ratioB+ratioBerr

oplot, density2, ratior, psym = 4, thick = 3, color = colors.black
errplot, density2,ratior-ratiorerr, ratior+ratiorerr
spearman = r_correlate(density2, ratior)
;xyouts, 50, 26, string(spearman(0)), charthick = 3

start=[-1E11,20]
result = MPFITFUN('linear',density2, ratior,ratiorerr, start,perror =test,/quiet)
;;print, test
;oplot, density2,result(0)*density2 + result(1), color = colors.black
;oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
;oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black
;
;;;;;
;plot, density2, colorgrad, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.8Mpc", $
;  ytitle = "color gradient? 0=blue, 1=flat, 2=red",xrange=[50,160],$
;   charthick = 3, yrange=[-0.5,2.5],xstyle=1,ystyle = 1,$
;  xthick = 3, ythick = 3;
;spearman = r_correlate(density2, colorgrad)
;xyouts, 50, 2.2, string(spearman(0)), charthick = 3
;
;start=[-1E11,1.5]
;result = MPFITFUN('linear',density2, colorgrad,err, start,perror =test,/quiet)
;;print, test
;oplot, density2,result(0)*density2 + result(1), color = colors.black
;oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
;oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black
;;;
;plot, density2, lumr, psym = 2, thick = 3, xtitle = "ngals/(Mpc^2) inside 0.8Mpc", $
;  ytitle = "ICL r-band luminosity",yrange = [1E11,11E11],$
;   charthick = 3, ystyle = 1,$
;  xthick = 3, ythick = 3
;errplot, density2, lumr-lumrerr, lumr + lumrerr

;coeff1 = ROBUST_LINEFIT(density2,lumr, yfit, sig, coeff_sig)
;oplot, density2,yfit, thick =3
;print, "sigma", sig
;oplot, density2,yfit-sig, linestyle=2
;oplot, density2,yfit+sig, linestyle=2;

;start=[-1E11,7E11]
;result = MPFITFUN('linear',density2, lumr,lumrerr, start,perror =test,/quiet)
;print, test
;oplot, density2,result(0)*density2 + result(1), color = colors.black
;oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
;oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black

;start=[4E11]
;result = MPFITFUN('noslope',density2, lumr,lumrerr, start,perror =test)
;print, test
;oplot, density2,findgen(n_elements(density2)) - findgen(n_elements(density2)) + result(0), color = colors.blue
;oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
;oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black

;spearman = r_correlate(density2, lumr)
;xyouts, 50, 1E12, string(spearman(0)), charthick = 3
;;;
plot, density2, lumrzke, psym = 2, thick = 3, xtitle = "Number of Galaxies", $
  ytitle = "ICL r-band Luminosity (solar luminosities)",yrange = [0,8E11],$
   charthick = 3,xstyle=1, ystyle = 1,xrange=[50,300],$
  xthick = 3, ythick = 3
errplot, density2, lumrzke-lumrzkeerr, lumrzke + lumrzkeerr

spearman = r_correlate(density2, lumrzke)
xyouts, 50, 7E11, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',density2, lumrzke,lumrzkeerr, start,perror =test,/quiet)
;print, test
oplot, density2,result(0)*density2 + result(1), color = colors.black
oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black

;;;
plot, density2, clusterrzke, psym = 2, thick = 3, xtitle = "Number of Galaxies", $
  ytitle = "Galaxy r-band luminosity (solar luminosities)",yrange = [0,5E12],$
   charthick = 3,xstyle=1, ystyle = 1,xrange=[50,300],$
  xthick = 3, ythick = 3
errplot, density2, clusterrzke-clusterrzkeerr, clusterrzke + clusterrzkeerr

spearman = r_correlate(density2, clusterrzke)
xyouts, 50, 4E12, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,7E11]
result = MPFITFUN('linear',density2, clusterrzke,clusterrzkeerr, start,perror =test,/quiet)
;print, test
oplot, density2,result(0)*density2 + result(1), color = colors.black
;oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
;oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black
;;;;
plot, density2, z, psym = 2, thick = 3, xtitle = "Number of Galaxies", $
  ytitle = "Redshift",yrange = [0,0.4],$
   charthick = 3,xstyle=1, ystyle = 1,xrange=[50,300],$
  xthick = 3, ythick = 3
;errplot, density2, clusterrzke-clusterrzkeerr, clusterrzke + clusterrzkeerr

spearman = r_correlate(density2, z)
xyouts, 50, 0.34, string(spearman(0), format='(F10.2)'), charthick = 3

start=[-1E11,0.3]
result = MPFITFUN('linear',density2, z, err, start,perror =test,/quiet)
;print, test
oplot, density2,result(0)*density2 + result(1), color = colors.black
oplot, density2, (result(0) -2*test(0))*density2 + (result(1) - 2*test(1)), color = colors.black
oplot, density2, (result(0)+2*test(0) )*density2 + (result(1) + 2*test(1)), color = colors.black


;---
;-------------------------------------------------------
;device, /close
;set_plot, mydevice

ps_close, /noprint, /noid


END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;cdcolor = 0.
;nocdcolor = 0.
;cdlum = 0.
;nocdlum = 0.
;cdcolorgrad = 0.
;nocdcolorgrad = 0.
;cdfracV = 0.
;nocdfracV = 0.
;cdfrac3V = 0.
;nocdfrac3V = 0.
;cdfracr = 0.
;nocdfracr = 0.
;cdfrac3r = 0.
;nocdfrac3r = 0.
;c1 = 0
;c2 = 0
;
;FOR j = 0, 9, 1 DO BEGIN
;   IF cd[j] EQ 1 THEN BEGIN
;       cdlum = cdlum + lumr[j] 
;       IF color[j] GE 0. THEN cdcolor = cdcolor + color[j] 
;       cdcolorgrad = cdcolorgrad + colorgrad[j] 
;       cdfracV = cdfracV + iclfracV[j]
;       cdfrac3V = cdfrac3V + ratioB[j]
;       cdfracr = cdfracr + iclfracr[j]
;       cdfrac3r = cdfrac3r + ratior[j]
;       c1 = c1 + 1
;   ENDIF ELSE BEGIN
;       nocdlum = nocdlum + lumr[j] 
;       nocdcolor = nocdcolor + color[j]
;       nocdcolorgrad = nocdcolorgrad + colorgrad[j] 
;       nocdfracV = nocdfracV + iclfracV[j]
;       nocdfrac3V = nocdfrac3V + ratioB[j]
;       nocdfracr = nocdfracr + iclfracr[j]
;       nocdfrac3r = nocdfrac3r + ratior[j]
;       c2 = c2 + 1
;   ENDELSE
;ENDFOR
;newx = [0,1]
;newcolory = [nocdcolor / c2,cdcolor / c1]
;newlumy = [nocdlum / c2,cdlum / c1]
;newcolorgrady = [nocdcolorgrad / c2,cdcolorgrad / c1]
;newfracVy = [nocdfracV / c2, cdfracV / c1]
;newfracry = [nocdfracr / c2, cdfracr / c1]
;newfrac3Vy = [nocdfrac3V / c2, cdfrac3V / c1]
;newfrac3ry = [nocdfrac3r / c2, cdfrac3r / c1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;cdcolor = 0
;nocdcolor = 0
;richcolor = 0
;norichcolor = 0
;rich3color = 0
;cdlum = 0
;nocdlum = 0
;richlum = 0
;norichlum = 0
;rich3lum = 0
;cdcolorgrad = 0
;nocdcolorgrad = 0
;richcolorgrad = 0
;norichcolorgrad = 0
;rich3colorgrad = 0
;cdfracV = 0
;nocdfracV = 0
;richfracV =0
;norichfracV =0
;rich3fracV =0
;cdfrac3V = 0
;nocdfrac3V = 0
;richfrac3V = 0
;norichfrac3V = 0
;rich3frac3V = 0
;cdfracr = 0
;nocdfracr = 0
;richfracr = 0
;norichfracr = 0
;rich3fracr = 0
;cdfrac3r = 0
;nocdfrac3r = 0
;richfrac3r = 0
;norichfrac3r = 0
;rich3frac3r = 0
;c1 = 0
;c2 = 0
;c3 = 0
;FOR j = 0, 9, 1 DO BEGIN
;   IF rich[j] EQ 1 THEN BEGIN
;       richlum = richlum + lumr[j] 
;       richcolor = richcolor + color[j] 
;       richcolorgrad = richcolorgrad + colorgrad[j] 
;       richfracV = richfracV + iclfracV[j]
;       richfrac3V = richfrac3V + ratioB[j]
;       richfracr = richfracr + iclfracr[j]
;       richfrac3r = richfrac3r + ratior[j]
;       c1 = c1 + 1
;   ENDIF 
;   IF rich[j] EQ 2 THEN begin
;       norichlum = norichlum + lumr[j] 
;       norichcolor = norichcolor + color[j]
;       norichcolorgrad = norichcolorgrad + colorgrad[j] 
;       norichfracV = norichfracV + iclfracV[j]
;       norichfrac3V = norichfrac3V + ratioB[j]
;       norichfracr = norichfracr + iclfracr[j]
;       norichfrac3r = norichfrac3r + ratior[j]
;       c2 = c2 + 1
;   ENDIF
;   IF rich[j] EQ 3 THEN begin
;       rich3lum = rich3lum + lumr[j] 
;       rich3color = rich3color + color[j]
;       rich3colorgrad = rich3colorgrad + colorgrad[j] 
;       rich3fracV = rich3fracV + iclfracV[j]
;       rich3frac3V = rich3frac3V + ratioB[j]
;       rich3fracr = rich3fracr + iclfracr[j]
;       rich3frac3r = rich3frac3r + ratior[j]
;       c3 = c3 + 1
;   ENDIF
;
;ENDFOR
;newrichx = [1,2,3]
;newrichcolory = [richcolor / c1,norichcolor / c2,rich3color / c3]
;newrichlumy = [richlum / c1, norichlum / c2,rich3lum / c3]
;newrichcolorgrady = [richcolorgrad / c1,norichcolorgrad / c2,rich3colorgrad / c3]
;newrichfracVy = [richfracV / c1,norichfracV / c2,rich3fracV / c3 ]
;newrichfracry = [richfracr / c1,norichfracr / c2,rich3fracr / c3 ]
;newrichfrac3Vy = [richfrac3V / c1,norichfrac3V / c2, rich3frac3V / c3]
;newrichfrac3ry = [richfrac3r / c1,norichfrac3r / c2, rich3frac3r / c3]


;average ellipticities of each clusters ICL profiles.
;and the implied change in the total ICL fluxes
;avrell=[0.56,0.63,0.40,0.60,0.61,0.30,0.57,0.,0.67,0.0]
;rellfactor = (((2.-avrell)^2.)/2.) / 2.
;rellfactor(6) = 0.58 ; make the added correction for the region which is added by hand.

;avBell=[0.53,0.61,0.53,0.51,0.60,0.28,0.60,0.,0.68,0.0]
;Bellfactor = (((2.-avBell)^2.)/2.) / 2.
;Bellfactor(6) = 0.61  ; make the added correction for the region which is added by hand.
