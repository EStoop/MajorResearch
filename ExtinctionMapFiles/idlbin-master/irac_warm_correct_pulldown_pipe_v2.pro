pro irac_warm_correct_pulldown_pipe_v2, dir_loc
;dir_loc = '/Users/jkrick/virgo/irac/s18p14/r35320064'

  cd, dir_loc
  
  for ch = 0, 1 do begin
     
     command_bcd =  "find . -name '" + strcompress("SPITZER_I" +string(ch + 1) + "*_bcd.fits' ",/remove_all) + " > "+strcompress("ch" +string(ch +1) + "_bcd.lis",/remove_all)
     spawn, command_bcd
     command_bunc =  "find . -name '" + strcompress("SPITZER_I" +string(ch + 1) + "*_bunc.fits' ",/remove_all) + " > "+strcompress("ch" +string(ch +1) + "_bunc.lis",/remove_all)
     spawn, command_bunc
     command_bimsk =   "find . -name '" + strcompress("SPITZER_I" +string(ch + 1) + "*_bimsk.fits' ",/remove_all) + " > "+strcompress("ch" +string(ch +1) + "_bimsk.lis",/remove_all)
     spawn, command_bimsk
     
     
     
     readcol, strcompress("ch" +string(ch +1) + "_bcd.lis",/remove_all), format='a', ifiles
     readcol, strcompress("ch" +string(ch +1) + "_bunc.lis",/remove_all), format='a', sfiles
     readcol, strcompress("ch" +string(ch +1) + "_bimsk.lis",/remove_all), format='a', mfiles

;
     
     pre='corr_'     
     ofiles=strarr(n_elements(ifiles))
     osfiles=strarr(n_elements(ifiles))
     omfiles=strarr(n_elements(ifiles))
     for j = 0L, n_elements(ifiles)-1 do begin
        
        pos1 = strpos(ifiles(j), 'S')
        ofiles(j) = strmid(ifiles(j),0, pos1 ) + pre + $
                    strmid(ifiles(j), pos1, 50)
        
                                ;
        
        pos2 = strpos(sfiles(j), 'S')
        osfiles(j) = strmid(sfiles(j),0, pos2 ) + pre + $
                     strmid(sfiles(j), pos2, 50)
        
                                ;
        
        pos3 = strpos(mfiles(j), 'S')
        omfiles(j) = strmid(mfiles(j),0, pos3 ) + pre + $
                     strmid(mfiles(j), pos3, 50)
        
     endfor
     
     openw,1,strmid(dir_loc,45) + '_corr_bcd.lis'
     printf,1,ofiles
     close,1
     
     openw,2,strmid(dir_loc,45)  + '_corr_bunc.lis'
     printf,2,osfiles
     close,2
     
     openw,3,strmid(dir_loc,45) + '_corr_bimsk.lis'
     printf,3,omfiles
     close,3
     
     
;; input files
     
     readcol,strmid(dir_loc,45) + '_corr_bcd.lis',ofiles,format='a'
     readcol,strmid(dir_loc,45) + '_corr_bunc.lis',osfiles,format='a'
     readcol,strmid(dir_loc,45) + '_corr_bimsk.lis',omfiles,format='a'
     
     
; running pulldown
     
     n = n_elements(ifiles)
     for i=0,n-1 do begin 
        
        print,'i =',i
        
        ifile = ifiles(i) & ofile = ofiles(i)
        sfile = sfiles(i) & osfile = osfiles(i)
        mfile = mfiles(i) & omfile = omfiles(i)
        ;print, ifile, sfile, mfile, ofile, osfile, omfile
        
        irac_warm_correct_pulldown_v2, ifile, sfile, mfile, ofile, osfile, $
                                       omfile, pthres, ft, backg, /usesig_weight, /verbose, /pthres
        
     endfor
     
  endfor ; for each channel

     
end
  
