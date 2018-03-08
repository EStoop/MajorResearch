pro catalog_rematch

restore, '/Users/jkrick/idlbin/objectnew.sav'
oldmatch = [82535,2939,34907,77318,78269,32192,78564,80705,86135,33084,80128,79564,33592,79726,78554,81081,2842,2324,81079,2758,79833,81673,80743,78720,2757,77617,80811,84634,80697,81502,82413,78987,77733,78979,82463,2753,86687,81344, 75878,75882,74994,75626,1975,74195,76204,983,75568,78009,78505,75246,74349,73072,77989,75012,2144,78429,78109,77891,77600,75808,74239,2745,73459,73566,74827,74589]
newmatch =[2034,2247,1258,1500,1087,1090,1507,1817,2035,2030,1145,2006,1244,1092,1088,1178,84553,15680,1821,15404,14740,2028,1176,1091,15688,1044,1177,1222,1168,2323,15647,1796,1055,2152,1235,32544,1250,1524,2314,2394,1983,1985,74320,12611,998,13342,2133,1060,1076,1759,921,879,1079,1986,13901,1505,2145,13913,1040,1496,1751,12793,2360,1474,1491,922]

print, n_elements(oldmatch), n_elements(newmatch)


objectnew[newmatch].mips24xcenter = objectnew[oldmatch].mips24xcenter
objectnew[newmatch].mips24ycenter = objectnew[oldmatch].mips24ycenter
objectnew[newmatch].mips24ra = objectnew[oldmatch].mips24ra
objectnew[newmatch].mips24dec = objectnew[oldmatch].mips24dec
objectnew[newmatch].mips24flux = objectnew[oldmatch].mips24flux
objectnew[newmatch].mips24fluxerr = objectnew[oldmatch].mips24fluxerr
objectnew[newmatch].mips24mag = objectnew[oldmatch].mips24mag
objectnew[newmatch].mips24magerr = objectnew[oldmatch].mips24magerr
objectnew[newmatch].mips24bckgrnd = objectnew[oldmatch].mips24bckgrnd
objectnew[newmatch].iter = objectnew[oldmatch].iter
objectnew[newmatch].sharp = objectnew[oldmatch].sharp
objectnew[newmatch].mips24match = objectnew[oldmatch].mips24match
objectnew[newmatch].mips24nndist = objectnew[oldmatch].mips24nndist
objectnew[newmatch].mips24matchdist = objectnew[oldmatch].mips24matchdist

objectnew[oldmatch].mips24xcenter = 0.0
objectnew[oldmatch].mips24ycenter = 0.0
objectnew[oldmatch].mips24ra = 0.0
objectnew[oldmatch].mips24dec = 0.0
objectnew[oldmatch].mips24flux = -99.0
objectnew[oldmatch].mips24fluxerr = -99.0
objectnew[oldmatch].mips24mag = 99.0
objectnew[oldmatch].mips24magerr = 99.0
objectnew[oldmatch].mips24bckgrnd = 0.0
objectnew[oldmatch].iter = 0.0
objectnew[oldmatch].sharp = 0.0
objectnew[oldmatch].mips24match = 0.0
objectnew[oldmatch].mips24nndist = 0.0
objectnew[oldmatch].mips24matchdist = 0.0


print, objectnew[newmatch[2]].mips24ra

save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'
end
