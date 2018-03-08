pro wheel_speed, planetname
  t1 = systime(1)

  restore, filename = '/Users/jkrick/irac_warm/irac_wheel_speeds/speeds.sav'

; FIELD2 = P-0718, 719, 720, 721
;FIELD4 = speeds
;FIELD5 = sclk
 
;------------------------------------------------------------------------
;  planetname = ['hatp2', 'HD20003','HD7924', 'HD45184', 'HD93385','XO3','HD209458', '55cnc']; ['hatp2', 'HD20003','HD7924', 'HD45184', 'HD93385', '55cnc'];
  planetinfo = create_planetinfo()
  chname = '2'
  apradius = 2.25

  for p = 0, n_elements(planetname) -1 do begin
     if chname eq '2' then aorname= planetinfo[planetname(p), 'aorname_ch2'] else aorname = planetinfo[planetname(p), 'aorname_ch1'] 
     basedir = planetinfo[planetname(p), 'basedir']
     
     dirname = strcompress(basedir + planetname(p) +'/') ;+'/hybrid_pmap_nn/')
     savefilename = strcompress(dirname + planetname(p) +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
     print, 'restoring ', savefilename
     restore, savefilename

     cenyrange = [14.5, 15.5]
     if planetname(p) eq '55cnc' then cenyrange = [13.0, 16.0]

     for a = 0, n_elements(aorname) - 1 do begin
        print, 'aorname', aorname(a)
       sclk_0 = planethash[aorname(a),'sclktime_0'] ; in seconds
        bmjdarr = planethash[aorname(a), 'bmjdarr']
                                ;turn bmjd into sclk
        elapsed_days= bmjdarr(n_elements(bmjdarr)-1) - bmjdarr(0)
        elapsed_secs = elapsed_days*24.*60.*60.
        sclk = dindgen(n_elements(bmjdarr))*elapsed_secs / n_elements(bmjdarr)
        sclk = sclk + sclk_0
        sclk_end = sclk(n_elements(sclk) -1 )
        
        
                                ;find the corresponding sclk in the wheel speed table
        sc = where(time gt sclk_0 and time lt sclk_end , sclkcount)
        print, sclkcount, time(sc(0)),time(sc(sclkcount - 1))
        
        if sclkcount gt 300 then begin  ; don't bother with the first 30 min AORs
           pw = plot(time(sc), rw1(sc),  color = 'silver', xtitle ='Sclk_time', ytitle = 'Wheel Speeds RPM', $
                     title = planetname(p) + ' ' + aorname(a), axis_style = 1 ,MARGIN = [0.15, 0.15, 0.20, 0.15]) 
           pw = plot(time(sc), rw2(sc), color = 'silver',/overplot, axis_style = 1 )
           pw = plot(time(sc), rw3(sc),   color = 'silver',/overplot, axis_style = 1 )
           pw = plot(time(sc), rw4(sc), color = 'silver',/overplot, axis_style = 1 )
        
           psx = plot(sclk, planethash[aorname(a),'xcen'], '1s',  axis_style = 0,sym_size = 0.2,   sym_filled = 1,  $
                      color = 'blue', /current,MARGIN = [0.15, 0.15, 0.20, 0.15], yrange = cenyrange)
           psy= plot(sclk, planethash[aorname(a),'ycen'], '1s', axis_style = 0, sym_size = 0.2,   sym_filled = 1,  color = 'green', /current, yrange = psx.yrange,MARGIN = [0.15, 0.15, 0.20, 0.15])
        
        
           a_psx = axis('y', target = psx, LOCATION = [max(pw.xrange),0,0], textpos = 1, title = 'Centroids')
           xaxis = axis('x', target = pw, LOCATION = [0, max(pw.yrange)], showtext = 0, tickdir = 1)

        endif

     endfor
  endfor
   
  print, 'time check', systime(1) - t1
     
end
  
