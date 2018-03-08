pro print_corrflux, planetname, apradius, chname

  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  
  dirname = strcompress(basedir + planetname +'/')       
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150226.sav',/remove_all) ;
  restore, savefilename

  startaor = 1
  stopaor = 20
  for a = startaor,stopaor, 2 do begin
     outname = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_'+aorname(a)+'.txt',/remove_all)
     openw, outlun,outname,/GET_LUN
     time = (planethash[aorname(a),'bmjdarr']) 
     corrflux =planethash[aorname(a),'corrflux_d']
     corrfluxerr =planethash[aorname(a),'corrfluxerr']

     for i = 0, n_elements(time) - 1 do begin
        printf, outlun, time(i), corrflux(i)/ median(corrflux), ' ', corrfluxerr(i)/ median(corrflux), format = '(D, F10.7 ,A,  F10.8)'
     endfor
     close, outlun
     free_lun, outlun
  endfor


end

