pro readcat
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)

close, /all

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/palomar/lfc/mass.2.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

restore, '/Users/jkrick/idlbin/object.sav'
fits_read, "/Users/jkrick/palomar/LFC/coadd_r.fits", data, header
count = 0
n = 0
mass = fltarr(n_elements(object.mass) )
z = fltarr(n_elements(object.mass) )
age = fltarr(n_elements(object.mass))
metal = fltarr(n_elements(object.mass))

model = strarr(n_elements(object.mass))


for i =0, n_elements(object.rmaga) - 1 do begin

 ;  if object[i].irxcenter gt 265.1298 and  object[i].irxcenter lt 265.130 then begin
 ;     if object[i].irycenter GT 69.05 and  object[i].irycenter lt 69.06 then begin
 ;        adxy, header, object[i].ra, object[i].dec, x, y
 ;        print,object[i].ra, object[i].dec, object[i].rmaga, object[i].rmagerra, x, y
 ;        print, object[i]

;      endif
;   endif


if object[i].mass gt 0 then begin
 count = count + 1
   mass[n] = object[i].mass
   z[n] = object[i].zphot
   age[n] = object[i].massage
   model[n] = strmid(object[i].model,62,7)

   if strcmp( strmid(object[i].model, 54,2) ,'22',2) eq 1 then metal[n] = 0.0001
   if strcmp( strmid(object[i].model, 54,2) ,'32',2) eq 1 then metal[n] = 0.0004
   if strcmp( strmid(object[i].model, 54,2) ,'42',2) eq 1 then metal[n] = 0.004
   if strcmp( strmid(object[i].model, 54,2) ,'52',2) eq 1 then metal[n] = 0.008
   if strcmp( strmid(object[i].model, 54,2) ,'62',2) eq 1 then metal[n] = 0.02
   if strcmp( strmid(object[i].model, 54,2) ,'72',2) eq 1 then metal[n] = 0.05
   n = n + 1
endif



;if object[i].massprob gt 0.1 and i lt 9000 then begin
;   print, "object", i, object[i].massprob, object[i].mass, strmid(object[i].model,62,3)

;endif

endfor

;print, count, " galaxies with mass > 0 in mass fitting"
mass = mass[0:n-2]
z = z[0:n-2]
age  = age[0:n-2]
metal = metal[0:n-2]
model = model[0:n-2]
sfr = strarr(n-1)

print, strmid(model[22], 3,4)
newmass = mass[sort(mass)]

plothist, alog10(newmass), xhist, yhist,/noprint,bin = 0.1, thick = 3, xthick=3,ythick=3, xtitle = "log(mass)"


;plot, xhist, yhist,thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "mass"
;title = "redshift distribution ",$
;xstyle = 1, xrange=[0,2], yrange=[0,1000], ytitle = "Number"

;print, object.masschi

;plotimage, xrange=[x-100,x+100],yrange=[y-100,y+100], bytscl(data, min = -0.1, max = 1) ,$; bytscl(data, min=1,max=3),$
 ;/preserve_aspect, /noaxes, ncolors=10

mass1 = mass
mass2 = mass
mass3 = mass
mass4 = mass
z1 = z
z2 = z 
z3 = z
z4 = z

plot, z, alog10(mass), psym =3, thick = 3, xthick=3,ythick=3,$
      xtitle = "redshift", ytitle = "log(mass)", charthick = 3, xrange=[0,3], yrange=[5,15]

j = 0
h = 0
k = 0
l = 0

print, n_elements(mass)

for m = 0, n_elements(mass) -1 do begin
   if age[m] le 1 then begin
      mass1[j] = mass(m)
      z1[j] = z(m)
      j = j + 1
   endif
   if age[m] le 3 and age[m] gt 1 then begin
      mass2[h] = mass(m)
      z2[h] = z(m)
      h = h + 1
   endif
   if age[m] le 5 and age[m] gt 3 then begin
      mass3[k] = mass(m)
      z3[k] = z(m)
      k = k + 1
   endif
   if age[m] le 8 and age[m] gt 5 then begin
      mass4[k] = mass(m)
      z4[l] = z(m)
      l = l + 1
   endif


   if strmid(model[m],0,3) eq "ssp" then sfr[m] = 1
;   if strmid(model[m],0,3) eq "exp" then sfr[m] = 2
   if strcmp( strmid(model[m],3,4), '0.01',4) eq 1 then sfr[m] = 3 
   if strcmp( strmid(model[m],3,4), '0.03',4) eq 1 then sfr[m] = 4 
   if strcmp( strmid(model[m],3,4), '0.05',4) eq 1 then sfr[m] = 5 
   if strcmp( strmid(model[m],3,4), '0.07',4) eq 1 then sfr[m] = 6 
   if strcmp( strmid(model[m],3,4), '0.1*', 3) eq 1 then sfr[m] = 7
   if strcmp( strmid(model[m],3,4), '0.3*', 3) eq 1 then sfr[m] = 8
   if strcmp( strmid(model[m],3,4), '0.5*', 3) eq 1 then sfr[m] = 9
   if strcmp( strmid(model[m],3,4), '0.7*', 3) eq 1 then sfr[m] = 10
   if strcmp( strmid(model[m],3,3), '1.*', 2) eq 1 then sfr[m] = 11
   if strcmp( strmid(model[m],3,3), '3.*', 2) eq 1 then sfr[m] = 12
   if strcmp( strmid(model[m],3,3), '5.*', 2) eq 1 then sfr[m] = 13
   
 ;  if sfr[m] eq 0 then print, model

;   if strmid(model[m],3,4) eq string("0.1_") then sfr[m] = 5 
;   if strmid(model[m],3,4) eq string(0.3_) then sfr[m] = 6
;   if strmid(model[m],3,4) eq string(0.5_) then sfr[m] = 7 
;   if strmid(model[m],3,4) eq string(0.7_) then sfr[m] = 8 
;   if strmid(model[m],3,4) eq string(1._d) then sfr[m] = 9 
;   if strmid(model[m],3,4) eq string(3._d) then sfr[m] = 10 
;   if strmid(model[m],3,4) eq string(5._d) then sfr[m] = 11 

endfor
print, j,h,k,l

mass1 = mass1[0:j-1]
mass2 =mass2[0:h-1]
mass3 = mass3[0:k-1]
mass4 = mass4[0:l-1]
z1 = z1[0:j-1]
z2 = z2[0:h-1]
z3 = z3[0:k-1]
z4 = z4[0:l-1]



oplot, z1, alog10(mass1), thick =3, psym = 1, color = colors.black
oplot, z2, alog10(mass2), thick =3, psym = 1, color = colors.blue
oplot, z3, alog10(mass3), thick =3, psym = 1, color = colors.red
oplot, z4, alog10(mass4), thick =3, psym = 1, color = colors.yellow


plot, alog10(mass), metal, thick =3, psym = 2, /ylog, yrange = [1E-4, 1E-1], xrange=[8,15]
plot, alog10(mass), sfr, psym = 2
plothist, alog10(metal), xhist, yhist,/noprint,bin = 0.1, thick = 3, xthick=3,ythick=3, xtitle = "log(Z)"
plothist, sfr, xhist, yhist,/noprint,bin = 1, thick = 3, xthick=3,ythick=3, xtitle = "sfr"
plot, sfr, alog10(metal), psym = 2
device, /close
set_plot, mydevice

END

