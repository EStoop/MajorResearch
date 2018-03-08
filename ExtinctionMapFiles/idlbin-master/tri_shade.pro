; File: tri_shade.pro
; Author: Erik Brisson

pro triangulated_surf

n = 32
seed = 1l
x = randomu(seed, n)
y = randomu(seed, n)
z = exp(-3* ((x-0.5)^2 + (y-0.5)^2) )

help, x, y, z
triangulate, x, y, tr

window,0,retain=2, title='Triangulation of irregular data',xsize=500,ysize=500
plot, x, y, psym=1, title='Triangulation of irregular data'
for i=0, n_elements(tr)/3-1 do begin
    t = [tr(*,i),tr(0,i)]
    plots, x(t), y(t)
endfor

zmax = max(z)
meshed = trigrid(x,y,z,tr)
meshed_shades = byte(255.0*meshed)
print, "zmax, min(meshed), max(meshed) = ", zmax, min(meshed), max(meshed)
print, "min(meshed_shades), max(meshed_shades)", min(meshed_shades), max(meshed_shades)

;meshedtoobig = trigrid(x,y,z,tr,missing=zmax+1)
;diff = meshedtoobig - meshed
;print, min(diff), max(diff)


window,1,retain=2, title='Mesh plot of irregular data',xsize=500,ysize=500
surface, meshed, title='Mesh plot of irregular data'

window,2,retain=2, title='Surface plot with explicit shades',xsize=500,ysize=500
shade_surf, meshed, title='Surface plot of irregular data', xtitle='X',   ytitle='Y', ztitle='Z', charsize=2.0, shades=meshed_shades

window,3,retain=2, title='Surface plot of irregular data',xsize=500,ysize=500
shade_surf, meshed, title='Surface plot of irregular data', xtitle='X',   ytitle='Y', ztitle='Z', charsize=2.0

end
