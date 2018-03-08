pro virgo_find_cluster
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

;the region of the potential cluster
filename = '/Users/jkrick/Downloads/cluster_jkrick.csv'
readcol, filename, ra, dec, u, err_u, g, err_g, r ,err_r, i ,err_i,z ,err_z, photoz, photozerr, photozcc2, photozerrcc2,photozd1, photozerrd1

a = where( photozerrcc2 lt 0.2)
plothist, photozcc2(a), bin = 0.05, thick = 3

;need to look at some other regions to make sure this isn't
;just the standard redshift distribution.

filename = '/Users/jkrick/Downloads/off1_jkrick.csv'
readcol, filename, ra, dec, u, err_u, g, err_g, r ,err_r, i ,err_i,z ,err_z, photoz, photozerr, photozcc2, photozerrcc2,photozd1, photozerrd1

avgphotzarr = photozcc2
avgphotzerrarr = photozerrcc2

;---------------------------------------

filename = '/Users/jkrick/Downloads/off2_jkrick_0.csv'
readcol, filename, ra, dec, u, err_u, g, err_g, r ,err_r, i ,err_i,z ,err_z, photoz, photozerr, photozcc2, photozerrcc2,photozd1, photozerrd1
avgphotzarr = [avgphotzarr, photozcc2]
avgphotzerrarr = [avgphotzerrarr, photozerrcc2]

;---------------------------------------
filename = '/Users/jkrick/Downloads/off3_jkrick.csv'
readcol, filename, ra, dec, u, err_u, g, err_g, r ,err_r, i ,err_i,z ,err_z, photoz, photozerr, photozcc2, photozerrcc2,photozd1, photozerrd1
avgphotzarr = [avgphotzarr, photozcc2]
avgphotzerrarr = [avgphotzerrarr, photozerrcc2]


a = where( avgphotzerrarr lt 0.2)
plothist,  avgphotzarr(a), xhist, yhist,bin = 0.05 , /noplot; ,/overplot, color = greencolor, peak = 1
yhist = yhist / 3

oplot, xhist, yhist, color = greencolor, psym = 10
;---------------------------------------

end
