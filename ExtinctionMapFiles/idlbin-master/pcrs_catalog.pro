pro pcrs_catalog

;readcol, '~/external/irac_warm/pcrs_catalog/pcrsguidestarcatalog.txt', star_id, junk, junk, Validity, Q, posEr, pErWk,  vMag, rightAscensn,  declination, prpMtnRA, prpMtnDc,  parllx, magEr,  raErr, declEr, mKER, mKED, plxEr, dOjbE, bkgEr, bstEr, P, M, X


;a = where(Q eq 0)
;b = where(Q eq 1)
;plothist, posEr(b)/ 18.25 * 24., xhist, yhist,  bin=1.0, /noprint,/noplot
;ph = barplot(xhist, yhist,  xtitle = 'PCRS Catalog Position Error at 2015.25 (mas)', ytitle = 'Number', fill_color = 'sky blue' , title = titlename, name = 'grade B')

;plothist, posEr(a)/ 18.25 * 24., xhist, yhist,  bin=1.0, /noprint,/noplot
;ph1 = barplot(xhist, yhist, overplot = ph, fill_color = 'red', name = 'grade A')

;L = legend(target = [ph, ph1], position = [140,3000], /data, /auto_text_color)

;-----------

;hip = where(M eq 0)
;tyc = where(M gt 0)

;plothist, posEr(tyc)/ 18.25 * 24., xhist, yhist,  bin=1.0, /noprint,/noplot
;ph2 = barplot(xhist, yhist,  xtitle = 'PCRS Catalog Position Error at 2015.25 (mas)', ytitle = 'Number', fill_color = 'sky blue' , title = titlename, name = 'Tycho')

;plothist, posEr(hip)/ 18.25 * 24., xhist, yhist,  bin=1.0, /noprint,/noplot
;ph3 = barplot(xhist, yhist, overplot = ph2, fill_color = 'red', name = 'Hipparcos')

;L = legend(target = [ph2, ph3], position = [140,3000], /data, /auto_text_color)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
readcol, '~/external/irac_warm/pcrs_catalog/pubdb_pcrs.unl', star_id, Validity, Q, posEr, pErWk,  vMag, rightAscensn,  declination, prpMtnRA, prpMtnDc,  parllx, magEr,  raErr, declEr, mKER, mKED, plxEr, dOjbE, bkgEr, bstEr, P, M, L, epoch, x, y, z, spt_ind, CNTR, delimiter = '|',/debug

hip = where(M eq 0)
tyc = where(M gt 0)

plothist, posEr(tyc)/ 18.25 * 24., xhist, yhist,  bin=1.0, /noprint,/noplot
ph4 = barplot(xhist, yhist,  xtitle = 'PCRS Catalog Position Error at 2015.25 (mas)', ytitle = 'Number', fill_color = 'sky blue' , title = titlename, name = 'Tycho')

plothist, posEr(hip)/ 18.25 * 24., xhist, yhist,  bin=1.0, /noprint,/noplot
ph5 = barplot(xhist, yhist, overplot = ph4, fill_color = 'red', name = 'Hipparcos')

L = legend(target = [ph4, ph5], position = [140,3000], /data, /auto_text_color)

a = where(vMag gt 7.1 and vMag lt 7.2 and rightAscensn gt 170 and rightAscensn lt 174,namag)
print, 'namag', namag
if namag gt 0 then begin
   for i = 0, namag - 1 do begin
      print, rightAscensn(a(i)), declination(a(i)), vMag(a(i))
   endfor
endif


end
