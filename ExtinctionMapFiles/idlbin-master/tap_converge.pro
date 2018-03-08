pro tap_converge

;want to visually look at the output of TAP to see if it converged.
;http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf


;grab this readcol line from the top of the ascii file, or just change it to the correct directory
tap_readcol,'/Users/jkrick/irac_warm/pcrs_planets/hat22/TAPmcmc_20140211_1051/ascii_ch2_phot_TAP_MCMC.ascii',P,i,adr,rdr,tmid,mu1,mu2,ecc,omega,yint,slope,sigr,sigw,format='d,d,d,d,d,d,d,d,d,d,d,d,d'

; mixing traceplots
;want these to look even, with no long lingering at bad jumps
count = findgen(n_elements(P))

mixplot = plot(count, i, xtitle = 'iteration', ytitle = 'inclination')
;mixplot = plot(count, adr, xtitle = 'iteration', ytitle = 'a/R*')
;mixplot = plot(count, rdr, xtitle = 'iteration', ytitle = 'rp/R*')
;mixplot = plot(count, tmid, xtitle = 'iteration', ytitle = 'tmid')
;;mixplot = plot(count, ecc, xtitle = 'iteration', ytitle = 'Eccentricity')
;;mixplot = plot(count, omega, xtitle = 'iteration', ytitle = 'Omega')
;mixplot = plot(count, sigr, xtitle = 'iteration', ytitle = 'sigma red noise')
;mixplot = plot(count, sigw, xtitle = 'iteration', ytitle = 'sigma white noise')
;mixplot = plot(count, mu1, xtitle = 'iteration', ytitle = 'Linear limb darkening')
mixplot = plot(count, mu2, xtitle = 'iteration', ytitle = 'Quad limb darkening')

;make density plots
plothist, i, xhist, yhist,/noplot,/autobin;, bin = 1.
t = barplot(xhist, yhist, xtitle = 'inclination', ytitle = 'Number')
junk = plot_quartiles(i, yhist)


plothist, adr, xhist, yhist,/noplot,/autobin;, bin = 0.2
t = barplot(xhist, yhist, xtitle = 'a/R*', ytitle = 'Number')
junk = plot_quartiles(adr, yhist)

plothist, rdr, xhist, yhist,/noplot,/autobin;, bin = 0.01
t = barplot(xhist, yhist, xtitle = 'rp/R*', ytitle = 'Number')
junk = plot_quartiles(rdr, yhist)

plothist, tmid, xhist, yhist,/noplot,/autobin;, bin = 0.001
t = barplot(xhist, yhist, xtitle = 'tmid', ytitle = 'Number')
junk = plot_quartiles(tmid, yhist)

plothist, sigr, xhist, yhist,/noplot,/autobin;, bin = 0.0005
t = barplot(xhist, yhist, xtitle = 'sigma red noise', ytitle = 'Number')
junk = plot_quartiles(sigr, yhist)

plothist, sigw, xhist, yhist,/noplot,/autobin;, bin = 0.00005
t = barplot(xhist, yhist, xtitle = 'sigma white noise', ytitle = 'Number')
junk = plot_quartiles(sigw, yhist)

plothist, mu1, xhist, yhist,/noplot,/autobin;, bin = 0.1
t = barplot(xhist, yhist, xtitle = 'Linear LD', ytitle = 'Number')
junk = plot_quartiles(mu1, yhist)

plothist, mu2, xhist, yhist,/noplot,/autobin;, bin = 0.1
t = barplot(xhist, yhist, xtitle = 'Quad LD', ytitle = 'Number')
junk = plot_quartiles(mu2, yhist)

plothist, omega, xhist, yhist,/noplot,/autobin;, bin = 0.1
t = barplot(xhist, yhist, xtitle = 'Omega', ytitle = 'Number')
junk = plot_quartiles(omega, yhist)

plothist, ecc, xhist, yhist,/noplot,/autobin;, bin = 0.1
t = barplot(xhist, yhist, xtitle = 'Eccentricity', ytitle = 'Number')
junk = plot_quartiles(ecc, yhist)

end

function plot_quartiles, i, yhist
 ; print, 'yhist inside', yhist
  m = median(i)
;  print, 'm', m, 'n', n_elements(i)
  ymax =  max(yhist)
;  ymin = min(yhist)
  t = plot([m,m], [0, ymax], color = 'red',/overplot, thick = 3)
;now want the 15.9% and 84.1% values
;start with sorting the yhist values
  sortedi = i[Sort(i)] 
  fn = uint(.159*(n_elements(i)))
  low = sortedi(fn)
;  print, 'low', fn,  low
  e4 = uint(.841*(n_elements(i)))
  high = sortedi(e4)
;  print, 'hight', e4, high
  t = plot([low,low], [0, ymax], 'r--3',/overplot)
  t = plot([high,high], [0, ymax], 'r--3',/overplot)
  

  return, 0

end
