pro tap_combine_plot, chname, planetname
  dirname = '/Users/jkrick/irac_warm/pcrs_planets/'+planetname + '/'
  plot_dir = [dirname + 'TAPmcmc_nnnp_ch'+chname+'/',dirname + 'TAPmcmc_nn_ch'+chname+'/',dirname + 'TAPmcmc_pmap_ch'+chname+'/']
  colorarr = ['red', 'blue', 'grey']
  offset = [0,0.004, 0.008]
  offset_r = [0.987, 0.990, 0.993]
  for method = 0, n_elements(plot_dir) -1 do begin
     print, 'working on ', plot_dir(method)
     tap_readcol,plot_dir(method)+"ascii_phased_data.ascii",p,f,rf,t,m,r,format='(d,d,d,d,d,d)'
;     print, plot_dir(method)+"ascii_phased_model.ascii"
     tap_readcol,plot_dir(method)+"ascii_phased_model.ascii",t2,p2,mf2,format='(d,d,d)'


     model = interpol(mf2[sort(p2)],p2[sort(p2)],p[sort(p)])
     yr = [min(f)-(15d0*stdev(r)),max(f)+(3d0*stdev(r))+ max(offset)]

     p1 =  plot(p[sort(p)]*24d0,model+r[sort(p)] + offset(method), '1o',  sym_filled = 1, color=colorarr(method),yrange=yr, $
                xrange=mm(p*24d0),xthick=3,ythick=3,xtitle='Hours from Mid Transit',ytitle='Relative Flux',/overplot)
     p2 = plot(p[sort(p)]*24d0,model + offset(method),color=colorarr(method),thick=3,overplot = p1)
     p3 = plot(p[sort(p)]*24d0,r[sort(p)] + offset_r(method),color = colorarr(method),'1o', sym_filled = 1,overplot = p1)
     print, r[sort(p)] + offset_r(method)
     l = polyline([min(p1.xrange),max(p1.xrange)], [ mean(r)+ offset_r(method), mean(r)+ offset_r(method)],/data, color = colorarr(method), thick = 3);min(f)-(8d0*stdev(r))+offset_r(method), min(f)-(8d0*stdev(r))+offset_r(method)
  endfor


  p1.save, dirname + 'allfluxes_binned_ch'+chname+'_tap.eps'

end
