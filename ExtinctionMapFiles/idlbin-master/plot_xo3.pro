pro plot_xo3

  aorname = [' ',  '20150001', '20150003','20150005','20150007','20150009','20150011','20150013','20150015','20150017','20150019',' ']
  
  ;;--------------------------------
  ;;simulated data
  ;;--------------------------------
  ;;real pmap, binx10, chisq
  real_depth = [0.001517,0.00744,0.001837,0.001604,0.00216,0.001605,0.001943,0.001927,0.00165,0.00205]
  real_err = [0.0000845,0.000129,0.000142,0.0000924,0.000102,0.0000917,0.00009459,0.00009885,0.0000946,0.0001218]

  ;;simulated pmap, binx10, chisq
  simul_depth = [0.001663838,0.001829,0.0017329,0.001747,0.00187,0.001645,0.001842,0.0020144,0.00159,0.0021799]
  simul_err = [0.0000947,0.0000915,0.0001138,0.000092519,0.000105,0.000104,0.0000899,0.000111,0.0000967,0.00036]

  ;;simulated pmap, nobin, chisq
  simul_nobin_depth = [0.001673,0.001679,0.0016423,0.0015999,0.0016926,0.00165,0.0016419,0.0016677,0.00164336,0.0016]
  simul_nobin_err = [0.000106,0.000107,0.000111,0.000109,0.000114,0.000104,0.000106,0.0001059,0.00011,0.000108]

  ;;real pmap, nobin, tap in rp/rstar
  real_tap = [.0386, alog10(-1), .0355, .0397]
  real_tap_err = [.0036,alog10(-1),.0055,.0067]
  
  ;;simul pmap, nobin, tap  in rp/rstar
  simul_tap = [.0391, alog10(-1), .0351, 0.0402]
  simul_tap_err = [.0034, alog10(-1), .0047, .0043]

  ;;simul pmap, bin10, tap
  simul_tap_bin10 = [0.0402,alog10(-1), 0.0358,0.0411]
  simul_tap_bin10_err = [0.0028,alog10(-1),0.004,0.0034]


  ps1 = errorplot(findgen(10), real_depth, real_err, '1tu', sym_filled = 1, xtitle = 'AOR name', $
                  ytitle = 'eclipse depth', yrange = [0.001, 0.003], xrange = [-1,10], xtickname = aorname,$
                  xtext_orientation =30, color = 'orange', errorbar_color = 'orange', name = 'real pmap chisq bin10')
  ps11 = plot(findgen(10), fltarr(10) + median(real_depth), color = 'orange', overplot = ps1)


  ps2 = errorplot(findgen(10), simul_depth, simul_err, '1s', sym_filled = 1, color = 'red', $
                  errorbar_color = 'red', name = 'simul pmap chisq bin10', overplot = ps1)
  ps21 = plot(findgen(10), fltarr(10) + median(simul_depth), color = 'red', overplot = ps1)

  ps3 = errorplot(findgen(4), real_tap^2, real_tap_err^2, '1o', sym_filled =1, color = 'blue',$
                  errorbar_color = 'blue', name = 'real pmap tap nobin', overplot = ps1)
  ps31 = plot(findgen(4), fltarr(4) + median(real_tap^2), color = 'blue', overplot = ps1)

  ps4 = errorplot(findgen(4), simul_tap^2, simul_tap_err^2, '1S', sym_filled =1, color = 'cyan',$
                  errorbar_color = 'cyan', name='simul pmap tap nobin', overplot = ps1)
  ps41 = plot(findgen(4), fltarr(4) + median(simul_tap^2), color = 'cyan', overplot = ps1)


  ps5 = plot(findgen(10), fltarr(10) + .001875, name = 'input depth', overplot = ps1)


  ps6 = errorplot(findgen(10), simul_nobin_depth, simul_nobin_err, '1s', sym_filled = 1, color = 'pink', $
                  errorbar_color = 'pink', name = 'simul pmap chisq nobin', overplot = ps1)
  ps61 = plot(findgen(10), fltarr(10) + median(simul_nobin_depth), color = 'pink', overplot = ps1)

  
  ps7 = errorplot(findgen(4), simul_tap_bin10^2, simul_tap_bin10_err^2, '1S', sym_filled =1, color = 'sea green',$
                  errorbar_color = 'sea green', name='simul pmap tap bin10', overplot = ps1)
  ps71 = plot(findgen(4), fltarr(4) + median(simul_tap_bin10^2), color = 'sea green', overplot = ps1)

  ls = legend(target = [ps1, ps2, ps6, ps3, ps4, ps7, ps5], position = [3.5, .0030], /data, /auto_text_color)
  ;;--------------------------------
  ;;real data
  ;;--------------------------------

  aorname = [' ', 'r46471424', 'r46471168','r46470912', 'r46470656', '46470400','r46469632', 'r46469120','r46468608', 'r46468352', 'r46468096',' ']
  ;;pmap reduction, bin10, chisq
  chisq_bin10 = [0.00166,0.001903,0.001314,0.001207,0.001359,0.001533,0.00156,0.001328,0.00133,0.001321]
  chisq_bin10_err = [0.000099,0.00013,0.000104,0.000104,0.0001,0.000102,0.000136,0.0001235,0.0001239,0.000135]

  ;;pmap reduction, nobin, chisq
  chisq_nobin =[0.00146,0.001572,0.0014204,0.00147,0.001514,0.001483,0.001402,0.001455,0.001459,0.00126]
  chisq_nobin_err = [0.00011,0.0001217,0.0001209,0.000109,0.000108,0.000102,0.0001204,0.000104,0.000132,0.000149]

  ;;pmap reduction, bin10, tap
  tap_bin10 = [.0418,.0418,.0368, .0359, .0354,0.0399, .0377,.0358, .0338,.0332]
  tap_bin10_err = [.0026,.0037, .0024, .0026, .0023, .0032,.0039, .0046, .0061,.0046]
  
  ;;pmap reduction, nobin, tap
  ;;running now
  tap_nobin = [.0416, .0413,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.0375,.0306]
  tap_nobin_err = [.0031,.0024,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.0036,.0063]
  
  ;;bliss, nobin, chisq
  bliss_nobin_chisq = [.00183,.001723,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.00148,.00154]
  bliss_nobin_chisq_err = [9.18E-05,.0001317,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1), 9.29E-05, 9.88E-05]

  ;;bliss, nobin, tap
  bliss_nobin_tap = [0.0382,.0387,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.0401,.0419]
  bliss_nobin_tap_err = [.0019,.0021,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.0022,.0028]

 ;;need to add the values that Stevenson himself got
  stevenson = [.00185, .00167,alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.00150,.00152]
  stevenson_err = [.00009,.00009, alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),alog10(-1),.00008,.00008]

  
  p1 = errorplot(findgen(10), chisq_bin10, chisq_bin10_err, '1tu', sym_filled = 1, xtitle = 'AOR name', $
                  ytitle = 'eclipse depth', yrange = [0.001, 0.003], xrange = [-1,10], xtickname = aorname,$
                  xtext_orientation =30, color = 'red', errorbar_color = 'red', name = 'pmap chisq bin10')
  
  p11 = plot(findgen(10), fltarr(10) + median(chisq_bin10), color = 'red', overplot = p1)

  p2 = errorplot(findgen(10), chisq_nobin, chisq_nobin_err,  '1o', name = 'pmap chisq nobin',$
                 sym_filled = 1, color = 'pink', errorbar_color = 'pink', overplot = p1)
  p21 = plot(findgen(10), fltarr(10) + median(chisq_nobin), color = 'pink', overplot = p1)
  
  p3 = errorplot(findgen(10), tap_nobin^2, tap_nobin_err^2,  '1o', name = 'pmap tap nobin',$
                 sym_filled = 1, color = 'cyan', errorbar_color = 'cyan', overplot = p1)
  p31 = plot(findgen(10), fltarr(10) + median(tap_nobin^2), color = 'cyan', overplot = p1)

  p4 = errorplot(findgen(10), tap_bin10^2, tap_bin10_err^2,  '1o', name = 'pmap tap bin10',$
                 sym_filled = 1, color = 'blue', errorbar_color = 'blue', overplot = p1)
  p41 = plot(findgen(10), fltarr(10) + median(tap_bin10^2), color = 'blue', overplot = p1)

  p5 = errorplot(findgen(10), bliss_nobin_chisq, bliss_nobin_chisq_err,  '1S', name = 'bliss chisq nobin',$
                 sym_filled = 1, color = 'green', errorbar_color = 'green', overplot = p1)
  p51 = plot(findgen(10), fltarr(10) + median(bliss_nobin_chisq), color = 'green', overplot = p1)

  p6 = errorplot(findgen(10), bliss_nobin_tap^2, bliss_nobin_tap_err^2,  '1S', name = 'bliss tap nobin',$
                 sym_filled = 1, color = 'lime', errorbar_color = 'lime', overplot = p1)
  p61 = plot(findgen(10), fltarr(10) + median(bliss_nobin_tap^2), color = 'lime', overplot = p1)

   p7 = errorplot(findgen(10), stevenson, stevenson_err,  '1S', name = 'Stevenson fit',$
                 sym_filled = 1,  color = 'black', errorbar_color = 'black', overplot = p1)
   p71 = plot(findgen(10), fltarr(10) + median(stevenson), color = 'black', overplot = p1)

  ;;and add a mean value of everyone's reduction.
  pmean = plot(findgen(10), fltarr(10) + .0015445, color = 'black', name = 'average all users', overplot = p1)
  
  ls = legend(target = [p1, p2, p3, p4, p5, p6,p7, pmean], position = [3.5, .0030], /data, /auto_text_color)

end
