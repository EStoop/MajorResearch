pro initial_drift_phot
  t1 = systime(1)

;all subarray AORs with PCRS peakup from Ocotober 2014 - current day =
;                                                        April 2015
  aorname = [ 'r50678528',  'r50665984',  'r51824640',  'r51816192',  'r51836416',  'r51843584',  'r51819776',  'r51833088',  'r51832320',  'r51820544',  'r51844096',  'r51832576',  'r51837184',  'r51828480',  'r51840000',  'r51823360', 'r51840256',  'r51841536',  'r51842816',  'r51823872',  'r51831040',  'r51840512',  'r51834880',  'r51835136',  'r51820800',  'r51835904',  'r51827968',  'r51816704',  'r51834368',  'r51829504', 'r51816448',  'r51826176',  'r51822848',  'r51832064',  'r51838976',  'r51815936',   'r51835648',  'r51828224',  'r51831552',  'r51829248',  'r52364288', 'r52364544',  'r52364800',  'r52365312',  'r50651392',  'r50653440',  'r50653184',  'r50652928',  'r50654976',  'r50652672',  'r50655232',  'r50652416',  'r50678016',  'r50676480',   'r49711360',  'r51831296',  'r51818496',  'r51827200',  'r51842304',  'r51818752',  'r51830016',  'r51817728',  'r51816960',  'r47037184',  'r47029248',  'r51842560',  'r51839488',  'r51817216',  'r51821568', 'r51837440',  'r51830528',  'r47030784',  'r47053824',  'r51837952',  'r51826688',  'r51830784',  'r51843328',  'r51836672',  'r51833600',  'r51833344',  'r51838464',  'r51825152',  'r51825920',  'r51832832',  'r51838720', 'r51837696',  'r51823104',  'r51883264',  'r51884032',    'r51834624',  'r51826432',  'r51841024',  'r51842048',  'r51819520',  'r51824384',  'r51821056',  'r51836160',  'r51827456',  'r51843072',  'r53522688',  'r53522432',  'r53522176',  'r53521920',  'r53521152',  'r53521664',  'r53523200',  'r53522944',  'r54315776',  'r54315520',  'r50651136',  'r50656000',  'r50651648', 'r50655744',  'r50649856',  'r50656768',  'r50650624',  'r50656256',  'r50650112',  'r50653696',  'r50655488',  'r50651904'] ;'r51822080',  'r51834112',





  c = [ '1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',  '2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2']


  dirname = '/Users/jkrick/external/initial_drift/'
  planethash = hash()
  startfits = 0L
  for a = 0, n_elements(aorname) - 1 do begin ;   
     dir = strcompress(dirname +aorname(a)  ,/remove_all)
     print, 'working on ',dir
     if FILE_TEST(dir) gt 0 then begin
        CD, dir                 ; change directories to the correct AOR directory
        
        command  = strcompress( 'find ch'+c(a)+"/bcd -name 'SPITZER*_bcd.fits' > "+dir+'/ch'+c(a) +'/bcdlist.txt')
        spawn, command
        command2 =  strcompress('find ch'+c(a)+"/bcd -name '*bunc.fits' > "+dir+ '/ch'+c(a) +'/bunclist.txt')
        spawn, command2
        
        readcol,strcompress(dir +'/ch'+c(a) +'/bcdlist.txt'),fitsname, format = 'A', /silent
        readcol,strcompress(dir+'/ch'+c(a) +'/bunclist.txt'),buncname, format = 'A', /silent
        for i =startfits, n_elements(fitsname) - 1  do begin 
           header = headfits(fitsname(i)) ;
           sclk_obs= sxpar(header, 'SCLK_OBS')
           frametime = sxpar(header, 'FRAMTIME')
           bmjd_obs = sxpar(header, 'BMJD_OBS')
           ch = sxpar(header, 'CHNLNUM')
           ronoise = sxpar(header, 'RONOISE') ; these are zeros
           gain = sxpar(header, 'GAIN')
           fluxconv = sxpar(header, 'FLUXCONV')
           exptime = sxpar(header, 'EXPTIME')
           aintbeg = sxpar(header, 'AINTBEG')
           atimeend = sxpar(header, 'ATIMEEND')
           naxis = sxpar(header, 'NAXIS')
           framedly = sxpar(header, 'FRAMEDLY')
           pid = sxpar(header, 'PROGID')
           ra_ref = sxpar(header, 'RA_REF')
           dec_ref = sxpar(header, 'DEC_REF')
           campaign = sxpar(header, 'CAMPAIGN')
           
           if c(a) eq '2' and frametime eq 2 then ronoise = 12.1
           if i eq startfits and naxis eq 3 then begin
              xarr = fltarr(64*(n_elements(fitsname)))
              yarr = xarr
              timearr =  dblarr(64*(n_elements(fitsname)))
              bmjd = timearr
           endif
           
           deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
           nt = dindgen(64)
           sclkarr = sclk_obs  + (deltatime*nt)/60./60./24.D ; 0.5*frametime + frametime*nt
           bmjdarr= bmjd_obs + (deltatime*nt)/60./60./24.D   ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
           
           ;;read in the files
           fits_read, fitsname(i), im, h
           fits_read, buncname(i), unc, hunc
           
           ;;run the centroiding 
           get_centroids_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                            x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                            x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                            xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                            xfwhm, yfwhm,  /WARM
           x_center = temporary(x3)
           y_center = temporary(y3)

           if naxis eq 3 then begin ; and i eq 0 then begin
              xarr[i*64] = x_center
              yarr[i*64] = y_center
              timearr[i*64] = sclkarr     
              bmjd[i*64] = bmjdarr
           endif
           
         endfor                  ; for each fits image
        
        keys =['ra', 'dec', 'pid', 'campaign', 'ch','xcen', 'ycen', 'timearr', 'aor', 'bmjdarr','exptime']
        values=list(ra_ref,  dec_ref, pid, campaign, ch, xarr, yarr,  timearr, aorname(a), bmjd,  exptime)
        planethash[aorname(a)] = HASH(keys, values)
        
     endif                      ; does that directory exits
  endfor                        ; for each AOR
  savename = dirname + 'centroids.sav'
  save, planethash, filename=savename
  print, 'saving planethash', savename
  print, 'time check', systime(1) - t1
  

end



;rm -rf 51843072
;rm -rf 53522688
;rm -rf 53522432
;rm -rf 53522176
;rm -rf 53521920
;rm -rf 53521152
;rm -rf 53521664
;rm -rf 53523200
;rm -rf 53522944
;rm -rf 54315776
;rm -rf 54315520
;rm -rf 50651136
;rm -rf 50656000
;rm -rf 50651648
;rm -rf 50655744
;rm -rf 50649856
;rm -rf 50656768
;rm -rf 50650624
;rm -rf 50656256
;rm -rf 50650112
;rm -rf 50653696
;rm -rf 50655488
;rm -rf 50651904
;
