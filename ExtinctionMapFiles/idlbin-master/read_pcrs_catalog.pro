pro read_pcrs_catalog, dirloc, star_id, Validity, Q, posEr, pErWk,  vMag, rightAscensn,  declination, prpMtnRA, prpMtnDc,  parllx, magEr,  $
                       raErr, declEr, mKER, mKED, plxEr, dObjE, bkgEr, bstEr, P, M, L, epoch, x, y, z, spt_ind, CNTR

                       
l = file_lines(dirloc)
a = strarr(l)
openr,unit,dirloc,/get_lun
readf,unit,a
;close,unit
free_lun, unit
asplit = strsplit(a[1:*],'|',/extract)
asplit2 = asplit.ToArray()
star_id = asplit2[*,0]
Validity = LONG(asplit2[*,1])
Q = LONG(asplit2[*,2])
posEr = DOUBLE(asplit2[*,3])
pErWk = DOUBLE(asplit2[*,4])
vMag = DOUBLE(asplit2[*,5])
rightAscensn = DOUBLE(asplit2[*,6])
declination = DOUBLE(asplit2[*,7])
prpMtnRA = DOUBLE(asplit2[*,8])
prpMtnDc = DOUBLE(asplit2[*,9])
parllx = DOUBLE(asplit2[*,10])
magEr = DOUBLE(asplit2[*,11])
raErr = DOUBLE(asplit2[*,12])
declEr = DOUBLE(asplit2[*,13])
mKER = DOUBLE(asplit2[*,14])
mKED = DOUBLE(asplit2[*,15])
plxEr = DOUBLE(asplit2[*,16])
dObjE = LONG(asplit2[*,17])
bkgEr = DOUBLE(asplit2[*,18])
bstEr = DOUBLE(asplit2[*,19])
P = LONG(asplit2[*,20])
M = LONG(asplit2[*,21])
L = LONG(asplit2[*,22])
epoch = FLOAT(asplit2[*,23])
X = DOUBLE(asplit2[*,24])
Y = DOUBLE(asplit2[*,25])
Z = DOUBLE(asplit2[*,26])
spt_ind = LONG64(asplit2[*,27])
CNTR = LONG64(asplit2[*,28])

RETURN
END
