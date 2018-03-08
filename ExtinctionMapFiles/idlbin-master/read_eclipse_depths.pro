pro read_eclipse_depths

  readcol, '/Users/jkrick/external/irac_warm/XO3/eclipse_depths_MCMC.txt', MCMCname, MCMCdepth, MCMCerrpos, MCMCerrneg, $
           format = '(A, F10.6, F10.6, F10.6)'
  err = fltarr(2,12)
  err[1,*] = MCMCerrpos
  err[0,*] = MCMCerrneg
  readcol, '/Users/jkrick/external/irac_warm/XO3/eclipse_depths_chisq.txt', chisqname, junk, junk, chisqdepth, chisqerr, $
           format = '(A, F10.6, F10.6, F10.6, F10.6)'
  print, 'stddev MCMC, chisq', stddev(MCMCdepth,/nan), stddev(chisqdepth,/nan)
  p = errorplot(findgen(n_elements(MCMCdepth)) + 1, MCMCdepth,err, '1rs', sym_filled =1, $
                errorbar_color = 'red', yrange = [.0007, .0023], xrange = [0, 13], name = 'MCMC'+string(1E6*stddev(MCMCdepth,/nan)),$
               xtitle = 'AOR', ytitle = 'Eclipse Depth', title = 'XO-3b Real')
  p1 = plot(findgen(15), fltarr(15) + 0.0015, overplot = p)
  p2 = errorplot(findgen(n_elements(MCMCdepth)) + 1, chisqdepth, chisqerr,'1bo', sym_filled = 1, $
                errorbar_color = 'blue',name = 'Chisq'+string(1E6*stddev(chisqdepth,/nan)), overplot = p)

  l = legend(target = [p,p2], position = [6, 0.0010], /data, /AUTO_TEXT_COLOR)
end
