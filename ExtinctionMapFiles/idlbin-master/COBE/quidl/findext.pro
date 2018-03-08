;=================================================================
;
pro findext,lng,lat,d,abs,amod
;pro findext,lng,lat,d,amod
;
;  procedure to retrieve the absorption in V from three-dimensional 
;  grids, based on the Galactic dust distribution of Drimmel & Spergel.
;
;	Input:	
;	lng	gal long (degrees)
;	lat	gal lat (degrees)
;	d	los distance (kpc)
;
;	Output:
;	abs	extinction in V band with rescaling
;	amod	extinction without rescaling
;	
; R. Drimmel	2002
;=================================================================
; load absorption grids and make rescaling maps
;---------------------------------------------------------
   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/avdisk.xdr"
   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/avdloc.xdr"
   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/avspir.xdr"
   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/avori.xdr"
   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/avori2.xdr"
;   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/rf_allsky.xdr"
   restore, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/rf_new.xdr"

; number od COBE pixels
      nsky=393216

; build skymaps of rescaling parameters for each component
      dfac=replicate(1.0,nsky) & sfac=replicate(1.0,nsky)
      lfac=replicate(1.0,nsky)
      index=where(ncomp eq 1)
      dfac(index)=rfac(index)
      index=where(ncomp eq 2)
      sfac(index)=rfac(index)
      index=where(ncomp eq 3)
      lfac(index)=rfac(index)


; define abs
num = n_elements(d)
abs = dblarr(num)
avloc = dblarr(num)
abspir = dblarr(num)
absdisk = dblarr(num)

; to radians
l = lng*!dtor
b = lat*!dtor

; Now for UIDL code:
; -find the index of the corresponding COBE pixel
; dimensions of sixpack = 768 x 512 (= 393216)
	vectoarr = lindgen(768,512)	;necessary because I reduced the 
					; arrays to vectors.
	incoor=fltarr(num,2)
	incoor(*,0)=lng(*)
	incoor(*,1)=lat(*)
	pxindex=coorconv(incoor,infmt='L',outfmt='P',inco='G',outco='R9')
	pix2xy,pxindex,x_out,y_out,res=9,/sixpack
	tblindex = vectoarr(x_out,y_out)
	
; Sun's coordinates (get from dprms)
xsun=-8.0
zsun=0.015

; calculate the maximum distance in the grid
dmax = replicate(100.,num)
if b ne 0. then dmax = .49999/abs(sin(b)) - zsun/sin(b)
if cos(l) ne 0. then dmax = dmax < (14.9999/abs(cos(l)) - xsun/cos(l))
if sin(l) ne 0. then dmax = dmax < 14.9999/abs(sin(l))

; replace distance with dmax when greater
r=d
index = where(d ge dmax,n)
if n ne 0 then r(index) = dmax(index)

; heliocentric cartesian coordinates
x = r*cos(b)*cos(l) 
y = r*cos(b)*sin(l)
z = r*sin(b) + zsun

; for stars in Solar neighborhood
i = where(abs(x) lt 1. and abs(y) lt 2., nloc, complement=j, ncomplement=nj)
if(nloc ne 0) then begin

; define the local grid
  dx=0.02 & dy=0.02 & dz=0.02
  nx = 101 & ny=201 & nz=51

; grid indices
  xi = x[i]/dx + float(nx - 1)/2.
  yj = y[i]/dy + float(ny - 1)/2.
  zk = z[i]/dz + float(nz - 1)/2.

; interpolate
  avloc[i] = interpolate(avori2,xi,yj,zk,missing=0.)
endif

k = where(abs(x) lt 0.75 and abs(y) lt 0.75, nloc, complement=m, ncomplement=nm)
if(nloc ne 0) then begin

; define the local grid
  dx=0.05 & dy=0.05 & dz=0.02
  nx = 31 & ny=31 & nz=51

; grid indices
  xi = x[k]/dx + float(nx - 1)/2.
  yj = y[k]/dy + float(ny - 1)/2.
  zk = z[k]/dz + float(nz - 1)/2.

; interpolate
  absdisk[k] = interpolate(avdloc,xi,yj,zk,missing=0.)
endif

; galacto-centric cartesian
x = x + xsun

; stars beyond local grid:

;larger orion arm grid 
if(nj ne 0) then begin

; calculate the allowed maximum distance for larger orion grid
  dmax = 100.
  if b ne 0. then dmax = .49999/abs(sin(b)) - zsun/sin(b)
  if cos(l) gt 0. then dmax = dmax < (2.374999/abs(cos(l)))
  if cos(l) lt 0. then dmax = dmax < (1.374999/abs(cos(l)))
  if sin(l) ne 0. then dmax = dmax < (3.749999/abs(sin(l)))

; replace distance with dmax when greater
  r1=d
  index = where(d ge dmax,n)
  if n ne 0 then r1(index) = dmax(index)

; galactocentric centric cartesian coordinates
  x1 = r1*cos(b)*cos(l) + xsun
  y1 = r1*cos(b)*sin(l)
  z1 = r1*sin(b) + zsun

; define the grid
  dx=0.05 & dy=0.05 & dz=0.02
  nx = 76 & ny=151 & nz=51

; grid indices
  xi = x1[j]/dx + 2.5*float(nx - 1)
  yj = y1[j]/dy + float(ny - 1)/2.
  zk = z1[j]/dz + float(nz - 1)/2.

; interpolate
  avloc[j] = interpolate(avori,xi,yj,zk,missing=0.)

endif

; define the grid
dx=0.2 & dy=0.2 & dz=0.02
nx = 151 & ny=151 & nz=51

; grid indices
xi = x/dx + float(nx - 1)/2.
yj = y/dy + float(ny - 1)/2.
zk = z/dz + float(nz - 1)/2.

; interpolate
abspir = interpolate(avspir,xi,yj,zk,missing=0.)
if(nm ne 0) then absdisk[m] = interpolate(avdisk,xi[m],yj[m],zk[m],missing=0.)

; apply rescaling factors
abs = dfac(tblindex)*absdisk + sfac(tblindex)*abspir + lfac(tblindex)*avloc

print,  dfac(tblindex),  sfac(tblindex), lfac(tblindex)

; without rescaling:
amod = absdisk + abspir + avloc

; done
return
end











