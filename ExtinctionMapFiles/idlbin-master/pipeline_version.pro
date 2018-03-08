pro pipeline_version

  dirname = '/Volumes/external/irac_warm/darks/superdarks/s2s/'
  cd, dirname
  ntotdarks = 0

  aorname = ['r38025728','r38025984','r38026240','r38026496','r38026752','r38055936','r38074880','r38095104','r38161152','r38169088','r38179328','r38393344','r38396928','r38409472','r38413056','r38690304','r38693888','r38709248','r38712832','r38811648','r38815232','r38928384','r38931968','r39106304','r39109888','r39123968','r39127552','r39144704','r39148288','r39206656','r39210240','r39220992','r39224576','r39340800','r39344384','r39354624','r39358208','r39427840','r39431424','r39453184','r39456768','r39535872','r39539456','r39553280','r39556864','r40258560','r40262144','r40278272','r40281856','r40322560','r40326144','r41019904','r41023488','r41084160','r41087744','r41179136','r41180416','r41232896','r41236480','r41254656','r41258240','r41406976','r41410560','r41434368','r41437952','r41605632','r41609216','r41623296','r41626880','r41680640','r41684224','r41862144','r41865728','r41888512','r41892096','r41995264','r41998848','r42013184','r42016768','r42040832','r42044416','r42069760','r42066176','r42084096','r42087680','r42157824','r42159104','r42554880','r42559744','r42601216','r42604800','r42772992','r42776576','r43954432','r43958016','r43994112','r44156672','r44153088','r44235520','r44239104','r44376064','r44379648','r44396800','r44400384','r44422400','r44425984','r44485120','r44488704','r44618240','r44621824','r44920576','r45048832','r45052672','r45178880','r45182464','r45212416','r45216000','r45250304','r45281536','r45376256','r45379840','r45419008','r45422592','r45523456','r45527040','r45562112','r45565952','r45649408','r45671424','r45696768','r45731328','r45752064','r45772544','r45805312','r45808640','r45827328','r45860608','r45863936','r45966336','r46014208','r46084864','r46102528','r46098944','r46500096','r46503680','r46926080','r46929664','r47076864','r47080448','r47584768','r47608064','r47611904','r48401408','r48426752','r48490752','r48566784','r48570368','r48582912','r48609024','r48620032','r48626432','r48648960','r48659712','r48752896','r48763648','r48801792','r48812544','r48828672','r48840192','r48855040','r48867072','r48872192','r48904192','r48921344','r48940288','r48953088','r48987392','r48998144','r49005056','r49015808','r49022208','r49036544','r49053440','r49064192','r49079040','r49089792','r49101056','r49112064','r49122816','r49133568','r49232128','r49241600','r50096896','r50107648','r49871872','r49882624','r50692096','r50702848','r50763264','r50774016','r50806016','r50824192','r50847744','r50886400','r50907648','r50918656','r50958336','r50975744','r51000320','r51021824','r51032576','r51047424','r51057920','r51074560','r51187456','r51198208','r51441920','r51452672','r51515392','r51526144','r51762432','r51782656','r51793408','r51803136','r51813888','r51848960','r51859712','r52000768','r52011520','r52039168','r52105984']

  pipeline = strarr(n_elements(aorname))
  date = fltarr(n_elements(aorname))
  bigim = fltarr(32, 32, 128, n_elements(aorname)) ; 
  
 for a = 0, n_elements(aorname) - 1 , 3 do begin

    command  = strcompress( 'find '+ aorname(a) + "/ch2/bcd/ -name 'SPITZER*0010*_bcd.fits' > 'bcdlist.txt'")
    spawn, command
    readcol,strcompress('bcdlist.txt'),fitsname, format = 'A', /silent

    header = headfits(fitsname) 
    pipeline(a)= sxpar(header, 'CREATOR')
    date(a) = sxpar(header, 'MJD_OBS')

    darkname = sxpar(header, 'SKDKRKEY')
    darkepid = sxpar(header, 'SDRKEPID')
    campaign = sxpar(header, 'CAMPAIGN')
    print, campaign, '  ', aorname(a)

    if file_test(strcompress('/Volumes/iracdata/flight/IWIC/S19repro/darks/' + campaign + '/cal/',/remove_all)) gt 0 then begin
       ;;compare this dark from the SHA to the S19.1 dark on iracdata

       fits_read, strcompress(aorname(a) + '/ch2/cal/SPITZER_I2_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), shadarkdata, shadarkheader

       ;;find the corresponding dark on iracdata

       fits_read, strcompress('/Volumes/iracdata/flight/IWIC/S19repro/darks/' + campaign + '/cal/00' + strmid(aorname(a),1) + '/IRAC.2.*skydark.fits',/remove_all), s19darkdata, s19darkheader
    

       diff = s19darkdata - shadarkdata
       fits_write, strcompress(dirname + '/diffdarks_' + aorname(a) + '.fits',/remove_all), diff, shadarkheader
    endif 

    ;;look at the histogram of pixel values before and after the pipeline
    ;;change
;    fits_read, strcompress(aorname(a) + '/ch2/cal/SPITZER_I2_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), darkdata, darkheader
;    darkdata1 = darkdata[*,*,1]
;    plothist, darkdata1, xhist, yhist, bin = 0.2,/nan, /noprint, /noplot
;    if pipeline(a) eq 'S18.18.0' then fc = 'red' else fc = 'cyan'
;    b = barplot(xhist, yhist, xtitle = 'MJy/sr', ytitle = 'Number', xrange = [10, 35], fill_color = fc, overplot = b)

 endfor

; pipearr = float(strmid(pipeline, 1, 4))
; MJD2DATE, date, year, month, day
; print, year(0), month(0), day(0)
; p = plot(date - date(0) + 23, pipearr, '1*', xtitle = 'Days since 2009-09-01', ytitle = "Pipeline Version in the Darks")

; p = plot([769 + 23,769 + 23], [18.0, 19.2], color = 'red', thick = 2, overplot = p)
; p = text(780 + 23, 18.6, "start good pmap", color = 'red', overplot = p, /data)

 
end
