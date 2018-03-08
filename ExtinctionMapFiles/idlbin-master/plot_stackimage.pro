pro plot_stackimage, planetname

  planetinfo = create_planetinfo()
  chname = planetinfo[planetname, 'chname']
  ra_ref = planetinfo[planetname, 'ra']
  dec_ref = planetinfo[planetname, 'dec']
  aorname = planetinfo[planetname, 'aorname']
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')

  z = pp_multiplot(multi_layout=[3,3], global_xtitle='Time (hrs) ',global_ytitle='Fractional Change in SB')


  for a = 0, 0 do begin
     ;savstack = strcompress(dirname + aorname(a) + '/bin_stack.sav')
     savstack = strcompress(dirname + '/bin_stack.sav')
     restore, savstack
     
     for x = 14, 16 do begin
        for y = 14, 16 do begin
           pix = bin_stack[x,y,*]
           normpix = pix / median(pix) - 1
 ;          xy=z.plot(snapshots[a].xcen, snapshots[a].ycen, '6r1.',color = 'black', thick = 2, xminor = 0, xtickinterval = 0.5)
           c = z.plot( findgen(n_elements(normpix)) * 0.00467062 , normpix, '1', yrange = [-0.3, 0.3]);, xrange = [8.5,9.5])  ;based on 2 0.1s subarray frames binned together
  ;         t = text(10, -.17, /data, strcompress('['+string(x-15)+',' + string(y-15) + ']' ))
        endfor  ;for y
     endfor ;for x

  endfor

end
