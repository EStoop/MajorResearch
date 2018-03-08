pro read_XO3csv


filename = '/Users/jkrick/external/irac_warm/XO3/XO3_r46466816_phot_ch2.txt'
XO3 = read_csv(filename, header = keys, count = count)

print, 'count', count
print, keys
help, /struct, XO3
print, 'printing', XO3.FIELD02

end
