pro grasil_color_color
!P.multi = [0,1,1]
;read in Eric's fluxes from Grasil elliptical models as a
;function of time.
ps_start, filename= '/Users/jkrick/Virgo/grasil_color_color_2.ps'

restore, '/Users/jkrick/Virgo/grasil/seddat.sav'
savefile = filepath('seddat.sav',root_dir=['/Users/jkrick/Virgo/grasil'])
sobj = obj_new('IDL_Savefile', savefile)
scontents = sobj->contents()
print, scontents.n_var
snames = sobj->names()
print, snames
print, dblst

bAB = bAB[1:*]
vAB = vAB[1:*]
i1AB = i1AB[1:*]
i2AB = i2AB[1:*]

plot, bAB - vAB, i1AB-i2AB, xrange = [-0.5, 1.5], xstyle = 1, yrange = [-1, 2], xtitle = 'B-V (AB)', ytitle = '[3.6-4.5](AB)', thick = 3, charthick=3, xthick = 3, ythick = 3;, title = 'Grasil models including dust'

name = [ '0.2', '0.4', '0.8', '1.5', '2.0', '3.0','4.0','5.0', '8.9','11','13']
;for i = 0, n_elements(bAB) -1 do begin
   xyouts,  bAB[0] - vAB[0], i1AB[0]-i2AB[0], name[0] + 'Gyr', charthick = 3, orientation = 30.0
   xyouts,  bAB[1] - vAB[1], i1AB[1]-i2AB[1], name[1] + 'Gyr', charthick = 3, orientation = 30.0
   xyouts,  bAB[2] - vAB[2], i1AB[2]-i2AB[2], name[2] + 'Gyr', charthick = 3, orientation = 30.0
   xyouts,  bAB[3] - vAB[3], i1AB[3]-i2AB[3], name[3] + 'Gyr', charthick = 3, orientation = 30.0
      xyouts,  bAB[6] - vAB[6], i1AB[6]-i2AB[6], name[6] + 'Gyr', charthick = 3, orientation = 30.0

   xyouts,  bAB[10] - vAB[10], i1AB[10]-i2AB[10], name[10] + 'Gyr', charthick = 3, orientation = 30.0
;endfor

;inner, outer
m87_B_V = [1.0, 0.825] - 0.12
m87_1_2 = [-0.56, -0.16]
err_B_V = [ 0.03, 0.075]
err_1_2 = [0.1 , 0.56]

;oplot, m87_B_V, m87_1_2, linestyle = 2
oploterror, m87_B_V, m87_1_2, err_B_V, err_1_2,psym = 3, thick = 3
xyouts, 0.8, -0.7, 'Inner', charthick = 3
xyouts, 0.74, -0.12, 'Outer', charthick = 3


ps_end
end
