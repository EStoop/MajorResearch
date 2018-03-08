pro superdark_seasonal
                                ;changing to read in individual dark
                                ;files to make one superdark for the
                                ;warm mission 
  expname = ['s0p1s','s0p4s', 's2s']
  chname = '2'

  for e = 2, 2 do begin         ; n_elements(expname) -1 do begin
;for now exptime can be 's0p1s' or 's2s'
     case expname(e) OF
        's0p1s': exptime_d = 0.1
        's0p4s':exptime_d = 0.4
        's2s': exptime_d = 2.0
     endcase
     
     
     
     dirname = '/Volumes/external/irac_warm/darks/superdarks/'
     expdir = dirname + expname(e)
     cd, expdir
     ntotdarks = 0
     if e eq 0 then aorname = ['r38712320','r38811136','r38814720','r38927872','r38931456','r39105792','r39109376','r39123456','r39127040','r39144192','r39147776','r39206144','r39209728','r39220480','r39224064','r39340288','r39343872','r39354112','r39357696','r39427328','r39430912','r39452672','r39456256','r39535360','r39538944','r39552768','r39556352','r40258048','r40261632','r40277760','r40281344','r40322048','r40325632','r41019392','r41022976','r41083648','r41087232','r41178624','r41180928','r41232384','r41235968','r41254144','r41257728','r41406464','r41410048','r41433856','r41437440','r41605120','r41608704','r41622784','r41626368','r41680128','r41683712','r41861632','r41865216','r41888000','r41891584','r41994752','r41998336','r42012672','r42016256','r42040320','r42043904','r42083584','r42087168','r42157312','r42159616','r42554368','r42559232','r42600704','r42604288','r42772480','r42776064','r43953920','r43957504','r43993600','r44152576','r44156160','r44235008','r44238592','r44375552','r44379136','r44396288','r44399872','r44421888','r44425472','r44484608','r44488192','r44617728','r44621312','r44920064','r45048320','r45052160','r45178368','r45181952','r45211904','r45215488','r45249792','r45253376','r45281024','r45284864','r45375744','r45379328','r45418496','r45422080','r45522944','r45526528','r45561600','r45565440','r45648896','r45652480','r45670912','r45682944','r45696256','r45699840','r45730816','r45734400','r45751552','r45755136','r45772032','r45775360','r45804800','r45808128','r45826816','r45830144','r45860096','r45863424','r45965824','r45969152','r46013696','r46017024','r46081024','r46084352','r46098432','r46102016','r46499584','r46503168','r46925568','r46929152','r47076352','r47079936','r47584256','r47587840','r47607552','r47611392','r48400896','r48404480','r48426240','r48429824','r48490240','r48493824','r48566272','r48569856','r48582400','r48593152','r48608512','r48619520','r48625920','r48637184','r48648448','r48659200','r48752384','r48763136','r48801280','r48812032','r48828160','r48839680','r48854528','r48866560','r48871680','r48882432','r48903680','r48920832','r48939776','r48953600','r48986880','r48997632','r49004544','r49015296','r49021696','r49036032','r49052928','r49063680','r49078528','r49089280','r49100544','r49111552','r49122304','r49133056','r49231616','r49241088','r49871360','r49882112','r50096384','r50107136','r50691584','r50702336','r50762752','r50773504','r50805504','r50816256','r50823680','r50834432','r50847232','r50857984','r50885888','r50896640','r50907136','r50918144','r50957824','r50976256','r50999808','r51010560','r51021312','r51032064','r51046912','r51057408','r51074048','r51084800','r51186944','r51197696','r51441408'] ;'r41930496','r41934080','r44660224',
     if e eq 1 then aorname = ['r38024704','r38024448','r38024960','r38025216','r38025472','r38055680','r38074624','r38094848','r38160896','r38168832','r38179072','r38393088','r38396672','r38409216','r38412800','r38690048','r38693632','r38708992','r38712576','r38811392','r38814976','r38928128','r38931712','r39106048','r39109632','r39123712','r39127296','r39148032','r39144448','r39206400','r39209984','r39220736','r39224320','r39340544','r39344128','r39354368','r39357952','r39427584','r39431168','r39452928','r39456512','r39535616','r39539200','r39553024','r39556608','r40258304','r40261888','r40278016','r40281600','r40322304','r40325888','r41019648','r41023232','r41083904','r41087488','r41178880','r41180672','r41232640','r41236224','r41254400','r41257984','r41406720','r41410304','r41434112','r41437696','r41605376','r41608960','r41623040','r41626624','r41680384','r41683968','r41861888','r41865472','r41888256','r41891840','r41995008','r41998592','r42012928','r42016512','r42040576','r42044160','r42083840','r42087424','r42157568','r42159360','r42554624','r42559488','r42600960','r42604544','r42772736','r42776320','r43954176','r43957760','r43993856','r44156416','r44152832','r44235264','r44238848','r44375808','r44379392','r44396544','r44400128','r44422144','r44425728','r44484864','r44488448','r44617984','r44621568','r44920320','r45048576','r45052416','r45178624','r45182208','r45212160','r45215744','r45250048','r45253632','r45281280','r45285120','r45376000','r45379584','r45418752','r45422336','r45523200','r45526784','r45561856','r45565696','r45649152','r45652736','r45671168','r45683200','r45696512','r45700096','r45731072','r45734656','r45751808','r45755392','r45772288','r45775616','r45805056','r45808384','r45830400','r45827072','r45860352','r45863680','r45966080','r45969408','r46013952','r46017280','r46084608','r46081280','r46102272','r46098688','r46499840','r46503424','r46925824','r46929408','r47076608','r47080192','r47588096','r47584512','r47607808','r47611648','r48401152','r48404736','r48426496','r48430080','r48490496','r48494080','r48570112','r48566528','r48582656','r48593408','r48608768','r48619776','r48626176','r48637440','r48648704','r48659456','r48752640','r48763392','r48801536','r48812288','r48828416','r48839936','r48854784','r48866816','r48882688','r48871936','r48903936','r48914688','r48940032','r48953344','r48987136','r48997888','r49015552','r49004800','r49021952','r49036288','r49053184','r49063936','r49078784','r49089536','r49100800','r49111808','r49122560','r49133312','r49231872','r49241344','r50096640','r50107392','r49871616','r49882368','r50691840','r50702592','r50763008','r50773760','r50805760','r50816512','r50823936','r50834688','r50847488','r50858240','r50886144','r50896896','r50907392','r50918400','r50958080','r50976512','r51000064','r51010816','r51021568','r51032320','r51047168','r51057664','r51074304','r51085056','r51187200','r51197952','r51441664','r51452416','r51515136','r51525888','r51769344','r51762176','r51782400','r51793152','r51802880','r51813632','r51848704','r51859456','r52000512','r52011264','r52038912','r52049664','r52105728','r52116480'] ;'r41930752','r41934336','r44660480',
     if e eq 2 then aorname = ['r38025728','r38025984','r38026240','r38026496','r38026752','r38055936','r38074880','r38095104','r38161152','r38169088','r38179328','r38393344','r38396928','r38409472','r38413056','r38690304','r38693888','r38709248','r38712832','r38811648','r38815232','r38928384','r38931968','r39106304','r39109888','r39123968','r39127552','r39144704','r39148288','r39206656','r39210240','r39220992','r39224576','r39340800','r39344384','r39354624','r39358208','r39427840','r39431424','r39453184','r39456768','r39535872','r39539456','r39553280','r39556864','r40258560','r40262144','r40278272','r40281856','r40322560','r40326144','r41019904','r41023488','r41084160','r41087744','r41179136','r41180416','r41232896','r41236480','r41254656','r41258240','r41406976','r41410560','r41434368','r41437952','r41605632','r41609216','r41623296','r41626880','r41680640','r41684224','r41862144','r41865728','r41888512','r41892096','r41995264','r41998848','r42013184','r42016768','r42040832','r42044416','r42069760','r42066176','r42084096','r42087680','r42157824','r42159104','r42554880','r42559744','r42601216','r42604800','r42772992','r42776576','r43954432','r43958016','r43994112','r44156672','r44153088','r44235520','r44239104','r44376064','r44379648','r44396800','r44400384','r44422400','r44425984','r44485120','r44488704','r44618240','r44621824','r44920576','r45048832','r45052672','r45178880','r45182464','r45212416','r45216000','r45250304','r45281536','r45376256','r45379840','r45419008','r45422592','r45523456','r45527040','r45562112','r45565952','r45649408','r45671424','r45696768','r45731328','r45752064','r45772544','r45805312','r45808640','r45827328','r45860608','r45863936','r45966336','r46014208','r46084864','r46102528','r46098944','r46500096','r46503680','r46926080','r46929664','r47076864','r47080448','r47584768','r47608064','r47611904','r48401408','r48426752','r48490752','r48566784','r48570368','r48582912','r48609024','r48620032','r48626432','r48648960','r48659712','r48752896','r48763648','r48801792','r48812544','r48828672','r48840192','r48855040','r48867072','r48872192','r48904192','r48921344','r48940288','r48953088','r48987392','r48998144','r49005056','r49015808','r49022208','r49036544','r49053440','r49064192','r49079040','r49089792','r49101056','r49112064','r49122816','r49133568','r49232128','r49241600','r50096896','r50107648','r49871872','r49882624','r50692096','r50702848','r50763264','r50774016','r50806016','r50824192','r50847744','r50886400','r50907648','r50918656','r50958336','r50975744','r51000320','r51021824','r51032576','r51047424','r51057920','r51074560','r51187456','r51198208','r51441920','r51452672','r51515392','r51526144','r51762432','r51782656','r51793408','r51803136','r51813888','r51848960','r51859712','r52000768','r52011520','r52039168','r52105984'] ;'r41931008','r41934592','r44660736',

;ugh 4 dimensions
     bigim_feb = fltarr(32, 32, 128, n_elements(aorname)*18)  ; XXXmaking this number up for now
     count_feb = 0
     bigim_may = fltarr(32, 32, 128, n_elements(aorname)*18)  ; XXXmaking this number up for now
     count_may = 0
     bigim_aug = fltarr(32, 32, 128, n_elements(aorname)*18)  ; XXXmaking this number up for now
     count_aug = 0
     bigim_Nov = fltarr(32, 32, 128, n_elements(aorname)*18)  ; XXXmaking this number up for now
     count_Nov = 0
     testpix = fltarr( n_elements(aorname)*18)


     for a = 0, n_elements(aorname) - 1 do begin
;        print, 'working on ', aorname(a)
        cd, expdir
       ;make a superskyflat for each subframe
        skyflatname = strcompress( aorname(a) + '/ch'+ chname + '/cal/*superskyflat*.fits',/remove_all)
;        print, 'skyflatname', skyflatname
        fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
        flat64 = fltarr(32,32,64)
        flatsingle = flatdata[*,*,0]
        for f = 0, 63 do flat64[*,*,f] = flatsingle


        cd, strcompress(aorname(a) + '/ch'+chname + '/bcd',/remove_all)
        command1 = 'ls *bcd.fits > bcdlist.txt'
        spawn, command1
        readcol,'bcdlist.txt',fitsname, format = 'A', /silent
        ntotdarks = ntotdarks + n_elements(fitsname)
        
 
;-----------------------
        for i = 0,  n_elements(fitsname) - 1 do begin

                                ;make sure I really pulled the correct
                                ;darks
           header = headfits(fitsname(i)) ;
           ch = sxpar(header, 'CHNLNUM')
           framtime = sxpar(header, 'FRAMTIME')
           naxis = sxpar(header, 'NAXIS')
           mjd_obs = sxpar(header, 'MJD_OBS')
           mjd2date, mjd_obs, year, month, day  ; figure out when the observations were taken

 ;          print, 'ch, naxis, exptime', ch, chname, naxis, framtime, exptime_d
              if chname eq ch and naxis eq 3 and framtime eq exptime_d then begin
                 fits_read, fitsname(i), data, header
        
                                ;back out the flux conversion
                 fluxconv = sxpar(header, 'FLUXCONV')
                 exptime = sxpar(header, 'EXPTIME')
                 data = data / fluxconv ; now in DN/s
                 data = data* exptime   ; now in DN
                                ; flip the image
                 data = reverse(data, 2)
                 
                                ;remove the flat
                 data = data * flat64
                 
                 darkname = sxpar(header, 'SKDKRKEY')
                 darkepid = sxpar(header, 'SDRKEPID')
                 framedelay = sxpar(header, 'FRAMEDLY')
                 aorkey = sxpar(header, 'AORKEY')
                 
                                ;remove the dark that was already used in the image
                 fits_read, strcompress('../cal/SPITZER_I'+chname+'_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), darkdata, darkheader
                 data = data + darkdata
                 
                                ;and put it together with the others
                                ;tomake the superdark


                 ;;ok instead of month now sort on year
                 CASE 1 of
                    ;;(month ge 1) and  (month le 3): begin
                    (year eq 2014): begin
                       bigim_feb(0, 0, 0, count_feb) = data 
                       count_feb = count_feb + 1
                    end
                     ;;(month ge 4) and  (month le 6): begin
                    (year eq 2011): begin
                       bigim_may(0, 0, 0, count_may) = data 
                       count_may = count_may + 1
                    end
                     ;;(month ge 7) and  (month le 9): begin
                    (year eq 2012): begin
                       bigim_aug(0, 0, 0, count_aug) = data 
                       count_aug = count_aug + 1
                    end
                     ;;(month ge 10) and  (month le 12): begin
                    (year eq 2013): begin
                       bigim_nov(0, 0, 0, count_nov) = data 
                       count_nov = count_nov + 1
                    end
                    else: asdf = 1
                 ENDCASE

;                 testpix[count] = data[11,11]
;                 count = count + 1
              endif             ; if the correct darks
           endfor               ; for each fitsname
        endfor                  ; for each AOR
;        print, 'final count', count, n_elements(aorname) * 18, ntotdarks
        bigim_feb= bigim_feb[*,*,*,0:count_feb-1]
        bigim_may= bigim_may[*,*,*,0:count_may-1]
        bigim_aug= bigim_aug[*,*,*,0:count_aug-1]
        bigim_nov= bigim_nov[*,*,*,0:count_nov-1]
        
        
;now make a median
        superdark_feb = median(bigim_feb, dimension = 4)
        superdark_may = median(bigim_may, dimension = 4)
        superdark_aug = median(bigim_aug, dimension = 4)
        superdark_nov = median(bigim_nov, dimension = 4)
        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_2014.fits',/remove_all), superdark_feb, header 
        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_2011.fits',/remove_all), superdark_may, header 
        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_2012.fits',/remove_all), superdark_aug, header 
        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_2013.fits',/remove_all), superdark_nov, header 
;        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_feb.fits',/remove_all), superdark_feb, header 
;        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_may.fits',/remove_all), superdark_may, header 
;        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_aug.fits',/remove_all), superdark_aug, header 
;        fits_write, strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_nov.fits',/remove_all), superdark_nov, header 

        ;;are these datasets drawn from the same distribution?
 

;make some difference images

        fits_read,  strcompress(dirname+ 'superdark_ch'+ chname +'_'+expname(e)+'_S19.fits',/remove_all), superdark, header 

        diff_feb = superdark - superdark_feb
        diff_may = superdark - superdark_may
        diff_aug = superdark - superdark_aug
        diff_nov = superdark - superdark_nov

        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_2014.fits',/remove_all), diff_feb, header 
        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_2011.fits',/remove_all), diff_may, header 
        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_2012.fits',/remove_all), diff_aug, header 
        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_2013.fits',/remove_all), diff_nov, header 
;        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_feb.fits',/remove_all), diff_feb, header 
;        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_may.fits',/remove_all), diff_may, header 
;        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_aug.fits',/remove_all), diff_aug, header 
;        fits_write, strcompress(dirname+ 'diffdark_ch'+ chname +'_'+expname(e)+'_nov.fits',/remove_all), diff_nov, header 


         
;------------------------------------------------

     endfor                     ; for each expname

end
  


;;and next somehow compare to the individual darks.
;;;darks are still in DN, so are quantized. => histograms won't help
;     fluxarr = fltarr(64*count)
;     c = 0
;     for i = 0, n_elements(bigim[0,0,0,*]) - 1 do begin
;                                ;divide superdark by individual darks
;;   divim = superdark / bigim[*,*,*,i]
;;   fits_write, strcompress(dirname+ 'divdark_'+expname(e)+'_'+ string(i) + '.fits',/remove_all), divim, header ; won't be the correct header
;        
;                                ;try getting rid of the flat and the zodi level
;        subim =    bigim[*,*,*,i] - superdark
;        subim = subim / flat64
; ;       subim = subim - mean(subim,/nan) ; get rid of the zodi (or any constant level in the image)
;        subim = subim*fluxconv
;        fits_write, strcompress(dirname+ 'subdark_'+expname(e)+'_'+ string(i) + '.fits',/remove_all), subim, header ; won't be the correct header
;        
;                                ;do aperture photometry at the center 
;        for j = 0, 63 do begin
;           aper, subim[*,*,j], 15.5, 15.5, flux, fluxerr, sky, skyerr, 1., 3., [3, 6],/nan, /exact,/flux,/silent
;           fluxarr(c) = flux
;           c = c + 1
;        endfor
;        
;     endfor
     
;;make a histogram of the empty fluxes
;     plothist, fluxarr, xhist, yhist, /autobin,/noplot;

;     start = [0.0, 2.5, 1000]
;     error = fltarr(n_elements(xhist)) + 0.1          ; uniform errors
;     result= MPFITFUN('mygauss',xhist,yhist, error, start) ;ICL
;;     h = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), thick = 2, color = 'green', /overplot)
;
;     ;normalize
;     yhist = yhist / float(n_elements(yhist));

;     if e eq 0 then h = plot(xhist, yhist, xtitle = 'Empty Aperture Flux(MJy/sr)', ytitle = 'Number', title = 'Test Super Darks', color = 'blue',; thick = 2, xrange = [-10, 10], name = '0.1s')
;     if e ne 0 then h2 = plot(xhist, yhist, color = 'green',thick = 2,/overplot, name = '2s')

;   l = legend(target = [h, h2], position = [9, 3.6], /data, /auto_text_color)


;aorname = [r38712320,r38811136,r38814720,r38927872,r38931456,r39105792,r39109376,r39123456,r39127040,r39144192,r39147776,r39206144,r39209728,r39220480,r39224064,r39340288,r39343872,r39354112,r39357696,r39427328,r39430912,r39452672,r39456256,r39535360,r39538944,r39552768,r39556352,r40258048,r40261632,r40277760,r40281344,r40322048,r40325632,r41019392,r41022976,r41083648,r41087232,r41178624,r41180928,r41232384,r41235968,r41254144,r41257728,r41406464,r41410048,r41433856,r41437440,r41605120,r41608704,r41622784,r41626368,r41680128,r41683712,r41861632,r41865216,r41888000,r41891584,r41930496,r41934080,r41994752,r41998336,r42012672,r42016256,r42040320,r42043904,r42083584,r42087168,r42157312,r42159616,r42554368,r42559232,r42600704,r42604288,r42772480,r42776064,r43953920,r43957504,r43993600,r44152576,r44156160,r44235008,r44238592,r44375552,r44379136,r44396288,r44399872,r44421888,r44425472,r44484608,r44488192,r44617728,r44621312,r44660224,r44920064,r45048320,r45052160,r45178368,r45181952,r45211904,r45215488,r45249792,r45253376,r45281024,r45284864,r45375744,r45379328,r45418496,r45422080,r45522944,r45526528,r45561600,r45565440,r45648896,r45652480,r45670912,r45682944,r45696256,r45699840,r45730816,r45734400,r45751552,r45755136,r45772032,r45775360,r45804800,r45808128,r45826816,r45830144,r45860096,r45863424,r45965824,r45969152,r46013696,r46017024,r46081024,r46084352,r46098432,r46102016,r46499584,r46503168,r46925568,r46929152,r47076352,r47079936,r47584256,r47587840,r47607552,r47611392,r48400896,r48404480,r48426240,r48429824,r48490240,r48493824,r48566272,r48569856,r48582400,r48593152,r48608512,r48619520,r48625920,r48637184,r48648448,r48659200,r48752384,r48763136,r48801280,r48812032,r48828160,r48839680,r48854528,r48866560,r48871680,r48882432,r48903680,r48920832,r48939776,r48953600,r48986880,r48997632,r49004544,r49015296,r49021696,r49036032,r49052928,r49063680,r49078528,r49089280,r49100544,r49111552,r49122304,r49133056,r49231616,r49241088,r49871360,r49882112,r50096384,r50107136,r50691584,r50702336,r50762752,r50773504,r50805504,r50816256,r50823680,r50834432,r50847232,r50857984,r50885888,r50896640,r50907136,r50918144,r50957824,r50976256,r50999808,r51010560,r51021312,r51032064,r51046912,r51057408,r51074048,r51084800,r51186944,r51197696,r51441408]
