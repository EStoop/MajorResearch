PRO ellipse
mydevice = !D.NAME
!p.multi = [0, 0, 1]
;SET_PLOT, 'ps'

a= 15.
b = 11.
xcenter = 0.
ycenter = 0.
PA = 20./!RADEG

npoints = 120
phi = 2 * !PI * (Findgen(npoints) / (npoints-1))
x = a*cos(phi)
y = b*sin(phi)

xprime = xcenter + (x*cos(PA)) - (y*sin(PA))
yprime = ycenter + (x*sin(PA)) + (y*cos(PA))

plot, xprime, yprime


xsortindex = sort(xprime)
xsort = xprime[xsortindex]
print, xsort
print, "then xrad = ",xsort(N_ELEMENTS(xprime) - 1)



end
