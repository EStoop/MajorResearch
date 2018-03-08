pro plot_test_bkgd
     savename = '/Users/jkrick/irac_warm/HD209458/testbkgd.sav'
     restore, savename

 ;------------------------------
     x = indgen(n_elements(fluxarr3_7),/long)
 ;    pf = plot(x, fluxarr3_7/ median(fluxarr3_7),'1s', sym_size = 0.2,   sym_filled = 1, color = 'black', $
;               xtitle = 'Frame Number', ytitle = 'Raw Flux', yrange = [0.96, 1.05])
 ;    pf = plot(x, fluxarr5_10/median(fluxarr5_10),'1s', sym_size = 0.2,   sym_filled = 1, color = 'red', overplot = pf)
 ;    pf = plot(x, fluxarr7_15/median(fluxarr7_15),'1s', sym_size = 0.2,   sym_filled = 1, color = 'blue', overplot = pf)
     ;-----

;     pb = plot(x, backarr3_7 ,'1s', sym_size = 0.2,   sym_filled = 1, color = 'black', $
;               xtitle = 'Frame Number', ytitle = 'Background', yrange = [0., .003])
;;     pb = plot(x, backarr5_10 + 0.0006, '1s', sym_size = 0.2,   sym_filled = 1, color = 'red', overplot = pb)
;     pb = plot(x, backarr7_15, '1s', sym_size = 0.2,   sym_filled = 1, color = 'blue', overplot = pb)

 ;------------------------------
;need some binning

     
     h = histogram(x, OMIN=om, binsize = 2*63L, reverse_indices = ri)
     bin_flux3_7 = dblarr(n_elements(h))
     bin_flux5_10 = bin_flux3_7
     bin_flux7_15 = bin_flux3_7
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
           
           meanclip, fluxarr3_7[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_flux3_7[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, fluxarr5_10[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_flux5_10[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, fluxarr7_15[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_flux7_15[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           c = c + 1
        endif
     endfor
     
     bin_flux3_7 = bin_flux3_7[0:c-1]
     bin_flux5_10 = bin_flux7_15[0:c-1]
     bin_flux7_15 = bin_flux7_15[0:c-1]
;------------------------------
     x = indgen(n_elements(bin_flux3_7),/long)
     pf = plot(x, bin_flux3_7/ median(bin_flux3_7),'1s', sym_size = 0.2,   sym_filled = 1, color = 'black', $
               xtitle = 'Frame Number', ytitle = 'Raw Flux', yrange = [0.98, 1.02])
 ;    pf = plot(x, bin_flux5_10/median(bin_flux5_10),'1s', sym_size = 0.2,   sym_filled = 1, color = 'red', overplot = pf)
     pf = plot(x, bin_flux7_15/median(bin_flux7_15),'1s', sym_size = 0.2,   sym_filled = 1, color = 'blue', overplot = pf)


end
