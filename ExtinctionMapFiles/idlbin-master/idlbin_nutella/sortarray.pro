pro sortarray

openw, outlun1, '/Users/jkrick/spitzer/mips/mips24/apexM1_step2/mosaic_extract_low.tbl',/get_lun
openw, outlun2, '/Users/jkrick/spitzer/mips/mips24/apexM1_step2/mosaic_extract_high.tbl',/get_lun

READCOL,'/Users/jkrick/spitzer/mips/mips24/apexM1_step2/mosaic_extract.tbl',F='I,F,F,F,F,I,I,F,F,F,F,A,F,F',num, ra, dec, x, y, fitx, fity, flux, snr, chi2,stat,deblend,aperture1,bad_pix1

sortindex = sort(flux)
;print, flux(sortindex)

for i = 0, (n_elements(num) - 1)/2, 1 do begin
printf, outlun1, format='(I10,F10.2,F10.2,F10.2,F10.2,I10,I10,F10.2,F10.2,F10.2,F10.2,A10,F10.2,F10.2)',num(sortindex(i)), ra(sortindex(i)), dec(sortindex(i)), x(sortindex(i)), y(sortindex(i)), fitx(sortindex(i)), fity(sortindex(i)), flux(sortindex(i)), snr(sortindex(i)), chi2(sortindex(i)),stat(sortindex(i)),deblend(sortindex(i)),aperture1(sortindex(i)),bad_pix1(sortindex(i))
endfor

for i = (n_elements(num) - 1)/2, n_elements(num) - 1, 1 do begin
printf, outlun2, format='(I10,F10.2,F10.2,F10.2,F10.2,I10,I10,F10.2,F10.2,F10.2,F10.2,A10,F10.2,F10.2)',num(sortindex(i)), ra(sortindex(i)), dec(sortindex(i)), x(sortindex(i)), y(sortindex(i)), fitx(sortindex(i)), fity(sortindex(i)), flux(sortindex(i)), snr(sortindex(i)), chi2(sortindex(i)),stat(sortindex(i)),deblend(sortindex(i)),aperture1(sortindex(i)),bad_pix1(sortindex(i))
endfor

close, outlun1
close, outlun2
free_lun, outlun1
free_lun, outlun2
end
