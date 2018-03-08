pro plume_size

araroi = [187.32094,187.27515,187.23817,187.31383,187.4071]
adecroi = [13.130524,13.084224,13.089371,12.87497,12.899801]

braroi = [187.13428,187.10259,187.27072,187.29272]
bdecroi = [13.079912,13.075608,12.843246,12.863825]

e1raroi = [187.40936,187.29532,187.31513,187.41925]
e1decroi = [12.880349,12.838504,12.761167,12.83201]

d1raroi = [186.88879,186.86213,186.85729,186.86701,186.88317,186.89125]
d1decroi =[13.22866,13.229458,13.190105,13.14446,13.136594,13.163355]

d2raroi=[186.88879,186.79582,186.80071,186.84598,186.85569,186.88317,186.86701,186.85729,186.86213]
d2decroi=[13.22866,13.254614,13.193228,13.168064,13.146031,13.136594,13.14446,13.190105,13.229458]



fits_read, '/Users/jkrick/Virgo/IRAC/s18p14/ch1/ch1_Combine-mosaic/mosaic.fits', data, header
naxis1 = fxpar(header, 'naxis1')
naxis2 = fxpar(header, 'naxis2')

adxy, header, e1raroi, e1decroi, px2, py2
      
photroi = Obj_New('IDLanROI', px2, py2) 
mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)
      
;apply the mask to make all data outside the roi = 0

onlyplume = data*(mask GT 0)

                                ;try turning 0's into  Nan's.
n = where(onlyplume eq 0)
                                ;print, 'n  ', n_elements(n)
nan = alog10(-1)
onlyplume(n) = nan

a = where(finite(onlyplume) gt 0 , goodarea)
a = where(finite(onlyplume)  lt 1, counta)
print, 'goodarea', goodarea     ;, '       nan', counta


end
