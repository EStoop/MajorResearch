command='rm lc.out '
spawn,command
command='f77 -O3 lightcurve.f -o lightcurve'
spawn,command
command='./lightcurve <lc.in > lc.out '
spawn,command
openr,1,'lc.out'
a=dblarr(2,2000)
readf,1,a
close,1
plot,a(0,*),a(1,*)-1.
