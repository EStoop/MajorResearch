function create_planetinfo

  planetinfo = hash()
  keys =['ra', 'dec','aorname_ch2','aorname_ch1','basedir','chname', 'utmjd_center', 'transit_duration', 'period', 'intended_phase', 'exptime','mask', 'maskreg','stareaor','plot_norm','plot_corrnorm']
chname = 0
;---
  planetname = 'WASP-80b'
  ra_ref = 303.167750
  dec_ref = -2.143944
  aorname_ch2 = ['r48360448', 'r48359936'] ;ch2
  aorname_ch1 = ['r48361472','r48360960']
  basedir = '/Users/jkrick/external/irac_warm/' 
  utmjd_center = 56124.91751D	;junk
  transit_duration= 127.8 ;min close but not perfect=
  period =   3.06785200  ;days 
  exptime = [0.4, 2.0]
  intended_phase = 0.0
  stareaor = 2
  plot_norm = 1.0
  plot_corrnorm = 1.0

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-38b'
  ra_ref = 243.95985208                    ; 16 15 50.3645 
  dec_ref = 10.03259056 ; +10 01 57.326 
  aorname_ch2 = ['r45676288', 'r45676032'] ;ch2
  aorname_ch1 = [' ']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56190.96147D	;55335.4205D
  transit_duration= 279.648 ;min 0.1942 days  taken from NSTED
  period =   6.871815  ;days 
  exptime = 0.4
  intended_phase = 0.5
  stareaor = 2
  plot_norm = 0.1039
  plot_corrnorm = 0.1025

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'fake_dec7'
  ra_ref = 36.712664
  dec_ref = 37.550367
  aorname_ch2 = ['r0000000000']   ;['r0000000001'] 
  aorname_ch1 = ['r0000000000']   ;['r0000000001'] 
  basedir = '/Users/jkrick/irac_warm/simulation_JI/' 
  utmjd_center = 55957.14809D
  transit_duration= 170.0; min  taken from exoplanet transit database
  period =    1.219867  ;days 
  exptime = 0.1
  intended_phase = 0.
  stareaor = 0
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-33b'
  ra_ref = 36.712664
  dec_ref = 37.550367
  aorname_ch2 = ['r45383424', 'r45384448', 'r45384704'] ;ch2
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55957.14809D    ;56029.12024D

  transit_duration= 170.0; min  taken from exoplanet transit database
  period =    1.219867  ;days 
  exptime = 0.4
  intended_phase = 0.
  stareaor = 3
  plot_norm = 0.175
  plot_corrnorm = 0.175

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'HD209458'
  ra_ref = 330.79488
  dec_ref = 18.8843175
;  aorname_ch2 = ['r48014336']
  
;  aorname_ch2 = ['r38703872','r38704128','r38704384','r38701312','r38703360','r38703616'] ;ch2 stares
  aorname_ch2 = ['r38703872','r38704128','r38704384','r38701312','r38703360','r38703616','r45188864','r45189120','r45189376','r45189632','r45189888','r45190144','r45190400','r45190656','r45190912','r45191168','r45191424','r45191680','r45191936','r45192192','r45192704','r45195264','r45192960','r45193216','r45193472','r45193984','r45193728','r45195520','r45194240','r45194496','r45194752','r45195008','r45196288','r45195776','r45197312','r45196032','r45196544','r45196800','r45197056','r45197568','r45197824','r45198080','r45192448']  ; all ch2
;  aorname_ch2 = ['r45188864','r45189120','r45189376','r45189632','r45189888','r45190144','r45190400','r45190656','r45190912','r45191168','r45191424','r45191680','r45191936','r45192192','r45192704','r45195264','r45192960','r45193216','r45193472','r45193984','r45193728','r45195520','r45194240','r45194496','r45194752','r45195008','r45196288','r45195776','r45197312','r45196032','r45196544','r45196800','r45197056','r45197568','r45197824','r45198080','r45192448'] ; just snaps
  aorname_ch1 = ['r50797056', 'r50800640', 'r50797312', 'r50797568', 'r50797824', 'r50798080', 'r50798336', 'r50798592', 'r50798848', 'r50799104', 'r50799360', 'r50799616', 'r50799872', 'r50800128', 'r50800384']  ; ch1 phase curve feb 2006, PID6XXX
  basedir = '/Users/jkrick/irac_warm/' 
  utmjd_center = 51818.05045D; 55938.54133D; 55942.06608D
  transit_duration=184.2 ; min
  period = 3.52474859           ;days intended_phase = 0.0
  exptime = 0.4
  intended_phase = 0
  stareaor = 5
  plot_norm = 0.46
  plot_corrnorm = 0.498

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
;---
  planetname = 'HD75289'
  ra_ref = 131.91819
  dec_ref = -41.73778
 aorname_ch2 = ['r52584192', 'r52616448', 'r52617728', 'r52619008', 'r52620288', 'r52621568',  'r52608256', 'r52616704', 'r52617984', 'r52619264', 'r52620544', 'r52621824', 'r52615680', 'r52616960', 'r52618240', 'r52619520', 'r52620800', 'r52622080', 'r52615936', 'r52617216', 'r52618496', 'r52619776', 'r52621056', 'r52622336', 'r52616192', 'r52617472', 'r52618752', 'r52620032', 'r52621312', 'r52622592']  ; all ch2  'r52655616', 'r52654080','r52652544','r52651008','r52649472','r52647936',
  basedir = '/Users/jkrick/external/irac_warm/' 
  utmjd_center = 57398.541D; 55938.54133D; 55942.06608D
  transit_duration=184.2 ; min
  period = 3.5091651          ;days intended_phase = 0.0
  exptime = 0.1
  intended_phase = 0
  stareaor = 0
  plot_norm = 1.0
  plot_corrnorm = 1.0

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-62b'
  ra_ref = 87.139641
  dec_ref =-63.988068
  aorname_ch2 = [ 'r49048576', 'r49048320', 'r49048064', 'r49047808', 'r49047552', 'r49047296', 'r49047040', 'r49046784', 'r49046528', 'r49046272', 'r49046016', 'r49045760', 'r49045504', 'r49045248', 'r49044992', 'r49044736', 'r49044480', 'r49044224', 'r49043968', 'r49043712', 'r49043456', 'r49043200', 'r49042944', 'r49042688', 'r49042432', 'r49042176', 'r49041920', 'r49041664', 'r49041408', 'r49041152', 'r49040896', 'r49040640', 'r49040384', 'r49040128', 'r49039872', 'r48695808']  ;'r48680448','r48708608', 'r48708352', 'r48708096', 'r48707840', 'r48707584', 'r48707328', 'r48707072', 'r48706816', 'r48706560', 'r48706304', 'r48702720', 'r48702464', 'r48702208', 'r48701952', 'r48701696', 'r48701440', 'r48701184', 'r48700928', 'r48700672', 'r48700416', 'r48697600', 'r48697344', 'r48697088', 'r48696832', 'r48696576', 'r48696320', 'r48696064', 'r48690688', 'r48690432', 'r48690176', 'r48689920', 'r48689664', 'r48689408', 'r48689152', 'r48688896', 'r48688640', 'r48687360', 'r48687104', 'r48686848', 'r48686592', 'r48686336', 'r48686080', 'r48685824', 'r48685568', 'r48685312', 'r48685056',
  aorname_ch1 = ['r48702976', 'r48680960']
;  aorname_ch2 = ['r48695808', 'r48680448']
;aorname = [  'r48680448'] ;['r48680960','r48702976']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56545.36259D
  transit_duration=228.7 ; min
  period = 4.411953          ;days intended_phase = 0.0
  exptime = 2.0
  intended_phase = 0.5
  stareaor = 82
  plot_norm = 0.0450
  plot_corrnorm = 0.0451

 values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---

  planetname = 'HD189733'
  ra_ref = 300.18185
  dec_ref = 22.709912
  aorname_ch1 = ['r41592320' ,'r41592832', 'r41591808','r41592576' ,'r41591552','r41592064', 'r41591296']  ;ch1
  aorname_ch2 = ['r38390016','r38390272','r38390528','r38390784','r38406656']  ;ch2 ;'r31756800','r36787968'
  basedir = '/Users/jkrick/irac_warm/' 
  utmjd_center =  55556.83649D; 55557.94753D 
  transit_duration=.07574*24.*60.
  period = 2.21857312
  intended_phase = 0
  exptime = 0.1
  stareaor = 5
  mask = 'yes'
  maskreg = fltarr(32,32)
  maskreg[13:16, 4:7] =  !Values.F_NAN
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---

  planetname = 'Kepler-5'
  ra_ref = 299.40663
  dec_ref = 44.034842
  aorname_ch1 = ['r38405632', 'r38406144' ]  ;ch1
  aorname_ch2 = ['r38404608', 'r38405120']  ;ch2 
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =  54955.40122
  transit_duration=0.
  period = 3.54846
  intended_phase = 0.5
  exptime = 12.
  stareaor = 0
  mask = 'no'
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
;---

  planetname = 'Kepler-17'
  ra_ref = 298.39583
  dec_ref = 47.815
  aorname_ch1 = ['r40251904', 'r40252416' ]  ;ch1
  aorname_ch2 = ['r40252160', 'r40251648']  ;ch2 
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   55185.17803
  transit_duration=0.
  period = 1.485710
  intended_phase = 0.0
  exptime = 12.
  stareaor = 0
  mask = 'no'
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'Kepler-9'
  ra_ref = 285.57398
  dec_ref = 38.401
  aorname_ch1 = ['r47052800', 'r47044608' ,'r47025920', 'r47044352','r47038976','r47043584','47053312','47043072']  ;ch1
  aorname_ch2 = ['r47027968', 'r47055104','r47043328', 'r47054848', 'r47056896', 'r47054592', 'r47028736', 'r47054336']  ;ch2 
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =  55072.93381
  transit_duration=0.
  period = 19.22418
  intended_phase = 0.0
  exptime = 2.
  stareaor = 0
  mask = 'no'
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'KOI069' ;kepler 93
  ra_ref = 291.41798
  dec_ref = 38.672208
  aorname_ch1 = []  ;ch1
  aorname_ch2 = ['r39438336', 'r39438592','r41009920','r41010432']  ;ch2 
  basedir = '/Users/jkrick/irac_warm/' 
  utmjd_center =  55072.93381 ;;junk
  transit_duration=0.
  period = 4.727 ;days
  intended_phase = 0.0
  exptime = [12.,12,2.,2.]
  stareaor = 0
  mask = 'no'
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'KOI245'  ;kepler 37 don't know if b, c, or d
  ra_ref = 284.05941
  dec_ref = 44.51844
  aorname_ch1 = []  ;ch1
  aorname_ch2 = ['r39420928','r41009664']  ;ch2 
  basedir = '/Users/jkrick/irac_warm/' 
  utmjd_center =  55072.93381 ;;junk
  transit_duration=0.
  period = 13.36 ;days
  intended_phase = 0.0
  exptime = [12.,12,2.,2.]
  stareaor = 0
  mask = 'no'
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'KOI246'  ;kepler 37 don't know if b, c, or d
  ra_ref = 291.0323
  dec_ref = 49.040266
  aorname_ch1 = []  ;ch1
  aorname_ch2 = ['r39371008']  ;ch2 
  basedir = '/Users/jkrick/external/irac_warm/' 
  exptime = [12.]
  stareaor = 0
  mask = 'no'
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
;---
  planetname = 'HD158460'   ; !!!This is a standard star, ie should have a flat light curve
  ra_ref = 261.42213
  dec_ref = 60.048479
  aorname_ch1 = ['r42051584', 'r42506496']  ;ch1
  aorname_ch2 = [ 'r45184256','r45184512','r45184768','r45185024','r45185280','r45185536','r45185792','r45186048','r45186304','r45186560','r45186816','r45187072','r45187328','r45187584','r45187840','r45188096','r45188352','r45188608'] ;0,1s snaps
;aorname_ch2 = ['r42051584', 'r42506496']; stares
;aorname_ch2 = ['r44499968','r44499712','r44497408','r44501504','r44501248','r44500992','r44500736', 'r44500480', 'r44500224','r44499456', 'r44499200','r44498944','r44498688','r44498432','r44498176','r44497920','r44497664','r44497152']  ; think these are saturated but not sure.

 basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =  55556.83649D; junk
  transit_duration=.07574*24.*60.  ; junk
  period = 1.4  ; junk
  intended_phase = 0  ; junk
  exptime = 0.1
  stareaor = 0
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'HD20003'
  ra_ref = 46.908914
  dec_ref = -72.321879
  aorname_ch1 = ['r0']  ;ch1
  aorname_ch2 = ['r48408064','r48408320']  
   basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  ;utmjd_center =  [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
  utmjd_center =  56537.696D
  transit_duration=0
  period = 11.85
  intended_phase = 0.5
  exptime = 0.4
  stareaor = 1
  plot_norm = .34708
  plot_corrnorm = .34497

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'HD1461'
  ra_ref = 4.6758985
  dec_ref = -8.0535363
  aorname_ch1 = ['r0']  ;ch1
  aorname_ch2 = ['r48815872']  ;ch2 
   basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =  56537.696D ; junk
  transit_duration=0 ;junk
  period = 2.5 ;junk
  intended_phase = 0.5
  exptime = 0.1
  stareaor = 1
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'hatp2'
  ra_ref = 245.1515
  dec_ref = 41.04789
  aorname_ch1 = ['r0']  ;ch1
  aorname_ch2 = ['r46477312']  ;ch2 ;'r31756800','r36787968'
   basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  ;utmjd_center =  [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
  utmjd_center =  0D
  transit_duration=0
  period = 0.
  intended_phase = 0.5
  exptime = 0.4
  stareaor = 1
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'HD45184'
  ra_ref = 96.18233
  dec_ref = -28.780888
  aorname_ch1 = ['r0']  ;ch1
  aorname_ch2 = ['r46917120']  ;ch2 ;'r31756800','r36787968'
   basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  ;utmjd_center =  [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
  utmjd_center =  55944.25889D  ; junk
  transit_duration=0
  period = 3.2 ;junk
  intended_phase = 0.
  exptime = 0.1
  stareaor = 1
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'HD93385'
  ra_ref = 161.56298
  dec_ref = -41.464344
  aorname_ch1 = ['r0']  ;ch1
  aorname_ch2 = ['r48407296']  ;ch2 ;'r31756800','r36787968'
   basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  ;utmjd_center =  [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
  utmjd_center =  55944.2D  ; junk
  transit_duration=0.5 ; junk
  period =3.2 ; junk
  intended_phase = 0.
  exptime = 0.1
  stareaor = 1
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = '55Cnc'
  ra_ref = 133.14756
  dec_ref = 28.330094
;aorname = [ 'r39524608', 'r42000384', 'r43981056', 'r43981312', 'r43981568', 'r43981824', 'r48069888', 'r48070144', 'r48070400', 'r48070656', 'r48070912', 'r48071168', 'r48071424', 'r48071680', 'r48071936', 'r48072192', 'r48072448', 'r48072704', 'r48072960', 'r48073216', 'r48073472', 'r48073728'] ; all ch2
 aorname_ch2 = ['r43981056', 'r43981312', 'r43981568', 'r43981824','r48069888','r48070144','r48070400','r48070656','r48070912','r48071168','r48071424','r48071680','r48071936','r48072192','r48072448','r48072704','r48072960','r48073216','r48073472','r48073728']  ; ch2
 aorname_ch1 = ['r0000000000']
 basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  ;utmjd_center =  [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
  utmjd_center =  55947.20505D
  transit_duration=105.7
  period = 0.73654 
  intended_phase = 0.5
  exptime = 0.02
  stareaor = 0
  plot_norm = 3.87
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'GJ3470'
  ra_ref = 119.77445833
  dec_ref = 15.39152778
  aorname_ch1 = ['r56751360' ,'r50648576', 'r50658048','r56751616' ] 
  aorname_ch2 = ['r0000000000']
  basedir = '/Users/jkrick/external/dither_preaor/' 
  utmjd_center = 56089.9769D
  transit_duration=96.300  ;XXX bullshit
  period = 3.33671
  intended_phase = 0.0
  exptime =  0.4
  stareaor = 0
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'gj436'
  ra_ref = 175.5493
  dec_ref = 26.7036
  aorname_ch2 = ['r38702592' ,'r38808064', 'r40848128','r42614016' ] 
  aorname_ch1 = ['r38807296']
  basedir = '/Users/jkrick/external/irac_warm/trending/' 
  utmjd_center = 55221.47190D
  transit_duration=96.300
  period = 2.64385
  intended_phase = 0.0
  exptime = [0.1, 0.4, 0.4, 0.4, 0.4]
  stareaor = 0
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'GJ1214'
  ra_ref = 258.83076
  dec_ref = 4.961675
  aorname_ch2 = [ 'r42045952', 'r42046208', 'r42046464', 'r42046720', 'r42046976', 'r42047232', 'r42047488', 'r42047744', 'r42048000', 'r42048256', 'r42048512', 'r42048768', 'r42049024', 'r42049280', 'r42049536', 'r42050048', 'r42050304', 'r42050560', 'r42050816', 'r42051072', 'r42051328'];, 'r42052096', 'r42052352', 'r42052864']   ;Deming program
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 2455391.65402D
  transit_duration=70. ;total junk 
  period = 1.58040482
  intended_phase = 0.0
  exptime = [2.0]
  stareaor = 50
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'HD7924b'
  ra_ref = 20.496279
  dec_ref = 76.710305
;  aorname_ch2= ['r44605184','r46981888']; 2012
  ;aorname = ['r44605184' ] ;2011
  aorname_ch2 = ['r44605184', 'r46981632','r46981888'] ; both 2011 and 2012
   aorname_ch1 = ['r0000000000']
   basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 54727.49 - 0.16 ;don't know where this additional change comes from
  utmjd_center = 55866.43 ;from Stephen
  ;utmjd_center = 55866.28554017 ; from the Knutson Nov 2011 observation
  transit_duration=77.73 ;guess at 0.01 period
  period = 5.3978
  intended_phase = 0.0
  exptime = 0.1
  stareaor = 3
  plot_norm = 1.46
  plot_corrnorm = 1.46

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-13b'
  ra_ref = 140.10268
  dec_ref = 33.882288
  aorname_ch2 = ['r45676544','r45675520']  ; ,; from my PID
                                ;aorname = ['r47047168', 'r47037952'] ;ch2 from Dessert PID this is suposedly the primary
  aorname_ch1 = ['r47036928', 'r47038208'] ;ch1 from Dessert PID this is NOT a  primary
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56086.49239D
  transit_duration=243.4
  period = 4.353011
  intended_phase = 0.5
  exptime = 0.4
  stareaor = 2
  plot_norm = 0.0383
  plot_corrnorm = 0.0278

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'wasp12'
  ra_ref = 97.6366
  dec_ref = 29.6723
  aorname_ch1 = ['r28578560'] ;ch1
  aorname_ch2 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 54773.146181D
  transit_duration=175.7
  period = 1.091423
  intended_phase = 0.5
  exptime = 12.
  stareaor = 0
  plot_norm = 0.0228
  plot_corrnorm = 1.0   ; not there

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'hat3'
  ra_ref = 206.0941
  dec_ref = 48.0287
  aorname_ch2= ['r31749376'] 
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55301.31066D
  transit_duration=124.6
  period = 2.899736
  intended_phase = 0.5
  exptime = 2.
  stareaor = 0
  plot_norm = .0274
  plot_corrnorm = .0274   ; not there

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'WASP-15b'
  ra_ref = 208.92783
  dec_ref = -32.159837
;  aorname = ['r45675264', 'r45675776'] 
  aorname_ch2 = [ 'r45675264', 'r45675776'] 
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56184.45438D
  transit_duration=222.900000 ; min from nsted
  period = 3.752066
  intended_phase = 0.5
  exptime = 2
  stareaor = 2
  plot_norm = 1.
  plot_corrnorm = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'WASP-16b'
  ra_ref = 214.68297
  dec_ref = -20.275605
  aorname_ch2 =['r45674240','r45674496'] ;['r48705536', 'r48692736'];
  aorname_ch1 = ['r48705536', 'r48692736']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56578.27412; 56179.09319
  transit_duration=115.2
  period = 3.118601
  intended_phase = 0.5
  exptime = 2
  stareaor = 2
  plot_norm = 1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'HAT-P-22'
  ra_ref = 155.68122
  dec_ref = 50.128716
  aorname_ch2 = ['r45675008','r45674752' ] 
  aorname_ch1 =['r49227008', 'r49227264']; ['r48704512','r48693504']     ;
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =  56081.30088 ;56486.04060;
  transit_duration=172.2
  period = 3.21222
  intended_phase = 0.5
  exptime = 0.4
  stareaor = 2
  plot_norm = 0.1146;1.
  plot_corrnorm = 0.1146

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'upsandb'
  ra_ref = 24.199343
  dec_ref = 41.405457
  aorname_ch2= ['r50986240',	'r50985728',	'r50985472',	'r50984960',	'r50984448',	'r50983936',	'r50983424',	'r50982656',	'r50982400',	'r50981632',	'r50981376',	'r50981120',	'r50980864',	'r50980608',	'r50980352',	'r50980096',	'r50982144',	'r50981888',	'r50979840',	'r50979584',	'r50989568',	'r50989312',	'r50989056',	'r50988288',	'r50988544',	'r50988800',	'r50988032',	'r50987776',	'r50987520',	'r50987264',	'r50987008',	'r50986752',	'r50986496']

  basedir = '~/external/irac_warm/' 
  utmjd_center = 0.   ; doesn't transit, so this isn't defined
  transit_duration = 0.; doesn't transit, so this isn't defined
  period = 4.61711
  intended_phase = 0.0
  exptime = .02
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = '2massJ15150083'
  ra_ref =228.74719
  dec_ref = 48.800812
  aorname_ch1= [ 'r53292288', 'r53286656','r53286144'] ;dither 'r54488576',

  basedir = '~/external/dither_preaor/' 
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = '2massJ08173001'
  ra_ref =124.37314
  dec_ref = -61.916068
  aorname_ch1= [ 'r53288704', 'r53283584','r53283328']; dither 'r54488320',

  basedir = '~/external/dither_preaor/' 
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = '2massJ05120636'
  ra_ref = 78.026482
  dec_ref = -29.831637
  aorname_ch1= [ 'r53282816', 'r53291008','r53290752'] ;dither'r54488064',

  basedir = '~/external/dither_preaor/' 
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-103b'
  ra_ref = 249.314867
  dec_ref = 7.183353
  aorname_ch2= [ 'r53510656', 'r53514240', 'r53513984', 'r53513472', 'r53512960'] ;dither 'r54487552',, 'r54487808'
  aorname_ch1 =[ 'r53510912', 'r53519104', 'r53518336', 'r53518080', 'r53519872'] ;dither 'r54487040',, 'r54487296'
  basedir = '~/external/dither_preaor/' 
  ;;these are incorrect
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 12.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'EPIC_203776696b'
  ra_ref = 246.61125000
  dec_ref = -24.97000000
  aorname_ch2= [ 'r57769472'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC_204890128b'
  ra_ref = 244.141815
  dec_ref = -20.400534
  aorname_ch2= [ 'r57767936'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC_205084841b'
  ra_ref = 248.399277
  dec_ref = -19.479185
  aorname_ch2= [ 'r57766144'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC_210558622b'
  ra_ref = 60.79283333
  dec_ref = 16.34744722
  aorname_ch2= [ 'r58854912'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC_210731500b'
  ra_ref = 66.23620833       
  dec_ref = 18.82741667 
  aorname_ch2= [ 'r58852864'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC_210968143b'
  ra_ref = 59.15241667       
  dec_ref = 22.47263056  
  aorname_ch2= [ 'r58854144'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC_212154564b'
  ra_ref = 133.641182
  dec_ref =  23.13289
  aorname_ch2= [ 'r59268352'] 

  basedir ='/Volumes/Backup2/jk/irac_warm/k2/' 
 
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
;---
  planetname = 'epic202083828'
  ra_ref = 94.206066
  dec_ref = 24.596486
  aorname_ch2= [ 'r54486528','r54341120', 'r54341376','r54486784'] ;dither

  basedir ='~/external/dither_preaor/' 
  ;;these are incorrect
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 12.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'EPIC204897050'
  ra_ref = 240.42087
  dec_ref = -20.368781
  ;;aorname_ch2= ['r63203584', 'r63204096','r63204608']
  aorname_ch2= [ 'r63204096']

  basedir = '~/external/irac_warm/stauffer/' 
  utmjd_center = double(57912.39952464)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

  values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
  
  
  ;;---
  planetname = 'EPIC204117263'
  ra_ref = 242.75894
  dec_ref = -23.597297
  ;;aorname_ch2= ['r63203584', 'r63204096','r63204608']
  aorname_ch2= [ 'r63203840']
  
  basedir = '~/external/irac_warm/stauffer/' 
  utmjd_center = double(57915.3073228)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
   planetinfo[planetname] = HASH(keys, values)


   
;---
  planetname = 'epic201912552'
  ra_ref = 172.560133
  dec_ref = 7.587843
  aorname_ch2= ['r56416512', 'r56416000']

  basedir = '~/external/irac_warm/' 
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 2
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'epic201637175'
  ra_ref = 169.482657
  dec_ref = 2.618976
  aorname_ch2= ['r56418304', 'r56419328']

  basedir = '~/external/irac_warm/' 
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 2
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'epic201367065b'
  ra_ref = 172.33545
  dec_ref = -1.45494
  aorname_ch2= ['r53522176','r53521920', 'r53521152', 'r53521664']

  basedir = '~/external/irac_warm/' 
  utmjd_center = double(56813.4189)
  transit_duration = 2.54 * 60. ;minutes
  period = 10.05403
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'epic201367065c'
  ra_ref = 172.33545
  dec_ref = -1.45494
  aorname_ch2= [ 'r53522944', 'r53523200']

  basedir = '~/external/irac_warm/' 
  utmjd_center = 2456812.2786 -2400000.5  
  transit_duration = 3.53 * 60. ;minutes
  period = 24.6454
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
;---
  planetname = 'epic201367065d'
  ra_ref = 172.33545
  dec_ref = -1.45494
  aorname_ch2= [ 'r53522432', 'r53522688']

  basedir = '~/external/irac_warm/' 
  utmjd_center = 2456826.2239 -2400000.5  
  transit_duration = 3.93 * 60. ;minutes
  period = 44.5619
  intended_phase = 0.0
  exptime = 2.0
  stareaor = 0
  plot_norm = 1.
  Plot = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'HAT-P-8'
  ra_ref = 343.0414
  dec_ref = 35.447403
  aorname_ch2= ['r36785664', 'r48695040', 'r48698112',  'r48699136',  'r48700160',  'r48704000',  'r48705280', 'r48684800',  'r48695296',  'r48698368',  'r48699392',  'r48703232',  'r48704256',  'r48705792', 'r48694528',  'r48695552',  'r48698624',  'r48699648',  'r48703488',  'r48704768',  'r48706048', 'r48694784',  'r48697856',  'r48698880',  'r48699904',  'r48703744',  'r48705024', 'r48694272', 'r48694016', 'r48693760', 'r48693248', 'r48692992', 'r48692480', 'r48692224', 'r48691968', 'r48684544', 'r48684288', 'r48684032', 'r48683520', 'r48683008', 'r48682496', 'r48681984', 'r48681728','r48691712', 'r48691456', 'r48691200', 'r48690944'] ;ch2 snaps with 2 AORS on front from Knutson program

  aorname_ch1 =  ['r36787456'] ; ch1 knutson eclipse['r36784384',
;  aorname = [ 'r36789504','r36785664'] ; ch2 knutson eclipse
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   55207.80851;55210.88489;
  transit_duration=228.5
  period = 3.076378
  intended_phase = 0.5
  exptime = 2.
  stareaor = 1
  plot_norm = 0.0448
  Plot = 0.0440

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'WASP-14b'
  ra_ref = 218.27649
  dec_ref = 21.894575
  aorname_ch2 = ['r45428992', 'r45428736', 'r45428480', 'r45428224', 'r45426688','r45840128','r45839616','r45839104','r45838592','r45847040','r45846784','r45846528','r45846272','r45846016','r45843200','r45842944','r45842688','r45842432','r45842176','r45841920','r45841664','r45841408','r45845504','r45840896','r45845760','r45845504','r45845248','r45844992','r45844736','r45844480','r45844224','r45843968','r45843712','r45843456','r45840640','r45840384','r45839872','r45839360','r45838848','r45838336','r48688384','r48688128','r48687872','r48687616','r48683776','r48683264','r48682752','r48682240','r48681472','r48681216','r48680704'] ;ch2 stares and time ordered snaps
;  aorname_ch2 = [ 'r45428992', 'r45428736', 'r45428480', 'r45428224', 'r45426688','r45838592', 'r45840128', 'r45841408', 'r45842176', 'r45842944', 'r45844480', 'r45845248', 'r45846016', 'r45846784', 'r45839104', 'r45840896', 'r45841664', 'r45842432', 'r45843200', 'r45844736', 'r45845504', 'r45846272', 'r45847040', 'r45839616', 'r45841152', 'r45841920','r45843968','r45843712','r45843456','r45840640','r45840384','r45839872','r45839360','r45838848','r45838336','r48688384','r48688128','r48687872','r48687616','r48683776','r48683264', 'r48682752','r48682240','r48681472','r48681216','r48680704','r45842688', 'r45844224', 'r45844992', 'r45845760', 'r45846528']  ; ch2 stares and snaps
;  aorname_ch2 = [ 'r45838592', 'r45840128', 'r45841408', 'r45842176', 'r45842944', 'r45844480', 'r45845248', 'r45846016', 'r45846784', 'r45839104', 'r45840896', 'r45841664', 'r45842432', 'r45843200', 'r45844736', 'r45845504', 'r45846272', 'r45847040', 'r45839616', 'r45841152', 'r45841920','r45843968','r45843712','r45843456','r45840640','r45840384','r45839872','r45839360','r45838848','r45838336','r48688384','r48688128','r48687872','r48687616','r48683776','r48683264', 'r48682752','r48682240','r48681472','r48681216','r48680704']  ; ch2 'allstares_ch2' ,
;aorname_ch2 = ['r45428992', 'r45428736', 'r45428480', 'r45428224', 'r45426688'] ;ch2 stares
;'r45843968','r45843712','r45843456','r45840640','r45840384','r45839872','r45839360','r45838848','r45838336'

  aorname_ch1 = [ 'allstares_ch1'];'r45426944', 'r45427200', 'r45427456','r45427712','r45427968']  ;ch1
;  aorname_ch1 = ['r31760384', 'allstares_ch1'];'r45426944', 'r45427200', 'r45427456','r45427712','r45427968']  ;ch1

  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56020.23972; 56172.81485; transit ;56187.39924  ;one of the secondary times
  transit_duration=183.6
  period = 2.243752
  intended_phase = 0.0
  exptime = 2
  stareaor = 4
  plot_norm =  0.0583 ; for 2.25 apradius
  plot_corrnorm =  0.0580; for 2.25 apradius
;  plot_norm =  0.0599  ; for 2.5 apradius
;  plot_corrnorm =  0.0595; for 2.5 apradius
;  plot_norm = .0611; 0.0599  ; for 3.0 apradius
;  plot_corrnorm = .0602; 0.0596; for 3.0 apradius
  mask = 'yes'
  maskreg = fltarr(32,32)
  maskreg[4:7, 12:14] =  !Values.F_NAN

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'WASP-76b'
  ra_ref = 26.63274583
  dec_ref = 2.70057361
  aorname_ch2 = [ 'r62823168','r62824448','r62824704','r62824960','r62823424']
  aorname_ch1 = [ 'r62822656','r62823680','r62823936','r62824192','r62822912']  ; ch1
  basedir = '/Users/jkrick/external/irac_warm/' 
  utmjd_center = 56107.250
  transit_duration=226.5
  period = 1.809886
  intended_phase = 0.5
  exptime = 2
  stareaor =0;4
  plot_norm = 1
  plot_corrnorm = 1
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'wasp7'
  ra_ref = 311.04294
  dec_ref = -39.225186
  aorname_ch2 = [ 'r62823168','r62824448','r62824704','r62824960','r62823424']
  aorname_ch1 = [ 'r62822656', 'r62823680','r62823936','r62824192','r62822912'] 
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55374.29647
  transit_duration=226.5
  period = 4.954658
  intended_phase = 0.5
  exptime = 2
  stareaor =0;4
  plot_norm = 1
  plot_corrnorm = 1
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'hatp17'
  ra_ref = 324.53611
  dec_ref = 30.488196
  aorname_ch2 = [ 'r42625280']
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55963.75327
  transit_duration=243.4
  period = 10.338
  intended_phase = 0.5
  exptime = 2
  stareaor =0;4
  plot_norm = 1
  plot_corrnorm = 1
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'hatp26'
  ra_ref = 213.15669
  dec_ref = 4.059612
  aorname_ch2 = [ 'r42624768']
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 56026.13620
  transit_duration=147.3
  period = 4.234516
  intended_phase = 0.5
  exptime = 2
  stareaor =0;4
  plot_norm = 1
  plot_corrnorm = 1
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'HD97658'
  ra_ref = 168.63792
  dec_ref = 25.710523
  aorname_ch2 = [ 'r42608128','r48823552']
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55981.5450
  transit_duration=147.3
  period = 4.234516
  intended_phase =0.0
  exptime = 0.1
  stareaor =0;4
  plot_norm = 0.01785
  plot_corrnorm = 0.01785
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-52b'   ;simu;ated
  ra_ref = 23.232338
  dec_ref = 8.7620157
  aorname_ch2 = [ 'r0000000000', 'r0000000001', 'r0000000002', 'r0000000003', 'r0000000004', 'r0000000005']
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55793.18143
  transit_duration=59. ;min
  period = 1.749
  intended_phase =0.0
  exptime = 2
  stareaor =5
  plot_norm = 0.01785
  plot_corrnorm = .0177
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'hip57274'
  ra_ref = 176.17053
  dec_ref = 30.95822
  aorname_ch2 = [ 'r44273920']
  aorname_ch1 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 54801.015
  transit_duration=184.8
  period = 8.135
  intended_phase =0.0
  exptime = 0.4
  stareaor =0;4
  plot_norm = 1
  plot_corrnorm = 1
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'corot2'
  ra_ref = 291.7771
  dec_ref = 1.3837
  aorname_ch1 = [ 'r31774976']
  aorname_ch2 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center = 55125.09100
  transit_duration=136.0
  period = 1.742994
  intended_phase =0.5
  exptime = 2.
  stareaor =0;4
  plot_norm = .0224
  plot_corrnorm = .0224
  
   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)


;---
  planetname = 'hat6'
  ra_ref = 354.77395
  dec_ref = 42.466135
  aorname_ch1 = [ 'r31758592']  ; ch1
  aorname_ch2 = ['r31751680']   ; ch2
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   55458.85371D; 55451.14774D
  transit_duration=210.4 
  period = 3.852985
  intended_phase = 0.5
  exptime = 2.
  stareaor = 0
  plot_norm =.050
  plot_corrnorm = .050

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---------
  planetname = 'xo4'
  ra_ref = 110.38827
  dec_ref = 58.268128
  aorname_ch1 = [ 'r31766784']
  aorname_ch2 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   55180.50972D
  transit_duration=100. ; junk
  period = 4.125083
  intended_phase = 0.5
  exptime = 2.
  stareaor = 0
  plot_norm =0.046
  plot_corrnorm = 0.046

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---------
  planetname = 'XO3'
  ra_ref = 65.469554
  dec_ref = 57.817067
  aorname_ch1 = [ 'r0']
  aorname_ch2 = ['r46467072', 'r46471424', 'r46467840', 'r46471168', 'r46470144', 'r46470912', 'r46467584', 'r46470656', 'r46469376','r46470400', 'r46466816', 'r46469632', 'r46468864', 'r46469120', 'r46469888', 'r46468608', 'r46467328', 'r46468352','r46471680', 'r46468096']
;  aorname_ch2 = [ 'r46466816', 'r46468096', 'r46469376', 'r46470656', 'r46484224', 'r46486784','r31618816', 'r46467072', 'r46468352', 'r46469632', 'r46470912', 'r46484992', 'r46487040','r39117312', 'r46467328', 'r46468608', 'r46469888', 'r46471168', 'r46485504', 'r46487296','r39117568', 'r46467584', 'r46468864', 'r46470144', 'r46471424', 'r46486016', 'r46487552','r39117824', 'r46467840', 'r46469120', 'r46470400', 'r46471680',  'r46486528','r24291328']
  basedir = '/Users/jkrick/external/irac_warm/' 
;  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   54024.897D ; junk
  transit_duration=100. ; junk
  period = 3.1915239 
  intended_phase = 0.5
  exptime = 2.
  stareaor = 20
  plot_norm =1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---------
  planetname = 'simul_XO3'
  ra_ref = 65.469554
  dec_ref = 57.817067
  aorname_ch1 = [ 'r0']
  aorname_ch2 = ['r20150000', 'r20150001', 'r20150002', 'r20150003', 'r20150004', 'r20150005', 'r20150006', 'r20150007', 'r20150008', 'r20150009', 'r20150010', 'r20150011', 'r20150012', 'r20150013', 'r20150014', 'r20150015', 'r20150016', 'r20150017', 'r20150018', 'r20150019']
  basedir = '/Users/jkrick/external/irac_warm/' 
  utmjd_center =   54024.897D ; junk
  transit_duration=100. ; junk
  period = 3.1915239 
  intended_phase = 0.5
  exptime = 2.
  stareaor = 20
  plot_norm =1.
  plot_corrnorm = 1.

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'wasp43'
  ra_ref = 154.90806
  dec_ref = -9.8065685
  aorname_ch2 = [ 'r42615040']  ;ch2
  aorname_ch1 = [ 'r42614272'] ;ch1
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   55772.00350D
  transit_duration=69.6
  period = 0.813475
  intended_phase = 0.5
  exptime = 2.
  stareaor = 0
  plot_norm =.0342
  plot_corrnorm = .0342

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'WASP-63b'
  ra_ref = 94.33750000
  dec_ref = -38.32333333
  aorname_ch2 = [ 'r51822080','r51834112']  ;ch2
  aorname_ch1 = [ 'r51820800','r51835904'] ;ch1
  basedir = '/Users/jkrick/external/irac_warm/' 
  utmjd_center =   55921.1527D
  transit_duration=100 ; junk
  period =4.37809
  intended_phase = 0.5
  exptime = 0.4
  stareaor = 0
  plot_norm =1.0
  plot_corrnorm = 1.0

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'tres2'
  ra_ref = 286.80833333
  dec_ref = 49.31638889
  aorname_ch2 = [ 'r46866688','r46857472']  ;ch2
  aorname_ch1 = [ 'r46859264','r46857728'] ;ch1
  basedir = '/Users/jkrick/external/irac_warm/' 
  utmjd_center = 53957.1355D
  transit_duration=100 ; junk
  period =2.4706133738
  intended_phase = 1.0
  exptime = 2.0
  stareaor = 0
  plot_norm =1.0
  plot_corrnorm = 1.0

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)


;---
  planetname = 'NoisePixTest2'
  ra_ref = 359.99931
  dec_ref = 0.000694
  aorname_ch1 = [ 'r0000000000'] ;ch1
  aorname_ch2 = ['r0000000000']
  basedir = '/Users/jkrick/irac_warm/simulation_JI/' 
  utmjd_center =   56701.677D
  transit_duration=69.6  ; junk
  period = 0.813475       ;junk
  intended_phase = 0.5
  exptime = 0.1
  stareaor = 0
  plot_norm =0.3852
  plot_corrnorm = 0.3852

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)
;---
  planetname = 'exptime0p02'
  ra_ref = 359.99931
  dec_ref = 0.000694
  aorname_ch1 = [ 'r0000000000'] ;ch1
  aorname_ch2 = [ 'r46039552', 'r46039808'];, 'r46042112', 'r46042368', 'r46044672', 'r46044928', 'r46047232', 'r46047488']
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/' 
  utmjd_center =   56701.677D ; junk
  transit_duration=69.6         ; junk
  period = 0.813475       ;junk
  intended_phase = 0.5
  exptime = 0.1
  stareaor = 0
  plot_norm =0.3852
  plot_corrnorm = 0.3852

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)

;---
  planetname = 'pmap_star_ch2'
  ra_ref = 269.72771
  dec_ref = 67.793367
;not sorted for subarray
 aorname_ch1 = ['r1']
aorname_ch2 = [ 'r41913600', 'r41913856', 'r41914112', 'r41914368', 'r41914624', 'r41914880', 'r41915136', 'r41915392', 'r41915648', 'r41915904', 'r41916160', 'r41916416', 'r41916672', 'r41916928', 'r41917184', 'r41917440', 'r41917696', 'r41917952', 'r41918208', 'r41918464', 'r42023424', 'r42023680', 'r42023936', 'r42024192', 'r42024448', 'r42024704', 'r42024960', 'r42025216', 'r42025472', 'r42025728', 'r42025984', 'r42026240', 'r42026496', 'r42026752', 'r42027008', 'r42027264', 'r42027520', 'r42027776', 'r42028032', 'r42028288','r42779904', 'r42780160', 'r42780416', 'r42780672', 'r42780928', 'r42781184', 'r42781440', 'r42781696', 'r42781952', 'r42782208', 'r42782464', 'r42782720', 'r42782976', 'r42783232', 'r42783488', 'r42783744', 'r42784256', 'r42784512', 'r42784768', 'r42794496', 'r42794752', 'r42795008', 'r42795264', 'r42795520', 'r42799616', 'r42799872', 'r42800128', 'r42800384', 'r42800640', 'r42800896', 'r42801152', 'r42801408', 'r42801664', 'r42801920', 'r42802176', 'r42802432', 'r42802688', 'r42802944', 'r42803200', 'r42803456', 'r42803712', 'r42803968', 'r42804224', 'r42804480', 'r44169472', 'r44169728', 'r44169984', 'r44170240', 'r44170496', 'r44170752', 'r44171008', 'r44171264', 'r44171520', 'r44171776', 'r44172032', 'r44172288', 'r44172544', 'r44172800', 'r44173056', 'r44173312', 'r44173568', 'r44173824', 'r44174080', 'r44174336', 'r44174592', 'r44174848', 'r44175104', 'r44175360', 'r44175616', 'r44175872', 'r44176128', 'r44176384', 'r44244480', 'r44244736', 'r44251648', 'r44251904', 'r44252160', 'r44254976', 'r44255232', 'r44255488', 'r44255744', 'r44256000', 'r44256256', 'r44256512', 'r44256768', 'r44257024', 'r44257280',  'r44407808', 'r44408064', 'r44408320', 'r44408576', 'r44408832', 'r44409088', 'r44409344', 'r44409600', 'r44409856', 'r44410112', 'r44451584', 'r44451840', 'r44452096', 'r44452352', 'r44452608', 'r44452864', 'r44453120', 'r44453376', 'r44453632', 'r44453888', 'r44454144', 'r44454400', 'r44454656', 'r44454912', 'r44455168', 'r44455424', 'r44455680', 'r44455936', 'r44456192', 'r44456448', 'r44456704', 'r44456960', 'r44457216', 'r44457472', 'r44457728', 'r44457984', 'r44458240', 'r44458496', 'r44458752', 'r44459008', 'r44459264', 'r44459520', 'r44459776', 'r44460032', 'r44460288', 'r44460544', 'r44460800', 'r44461056', 'r44461312', 'r44461568', 'r44461824', 'r44462080', 'r44462336', 'r44462592', 'r44462848', 'r44463104', 'r44463360', 'r44463616', 'r44463872', 'r44464128', 'r45217536', 'r45217792', 'r45255424', 'r45255680', 'r45286656', 'r45286912', 'r45381376', 'r45381632', 'r45423872', 'r45424128', 'r45528576', 'r45528832', 'r45567232', 'r45567488', 'r45655040', 'r45655296', 'r45673472', 'r45673728', 'r45702400', 'r45702656', 'r45736960', 'r45737216', 'r45757696', 'r45757952', 'r45777664', 'r45777920', 'r45810432', 'r45810688', 'r45832448', 'r45832704', 'r45865728', 'r45865984', 'r45971456', 'r45971712', 'r46019328', 'r46019584', 'r46035712', 'r46035968', 'r46036224', 'r46036480', 'r46036736', 'r46036992', 'r46037248', 'r46037504', 'r46037760', 'r46038016', 'r46038272', 'r46038528', 'r46038784', 'r46039040', 'r46039296', 'r50836736']

;aorname_ch2 = ['r46039040', 'r46039296', 'r46037504', 'r46037760', 'r46038016', 'r46038272', 'r46038528', 'r46038784']

;full array or only pipe0
; 'r42532864',  'r42533120', 'r42533376', 'r42533632', 'r42533888',  'r42534144', 'r42534400', 'r42534656', 'r42534912','r42535168', 'r42535424', 'r42535680', 'r42535936', 'r42536192', 'r42536448', 'r42536704', 'r42536960', 'r42537216', 'r42537472', 'r42537728', 'r44265216', 'r44265472', 'r44265728', 'r44265984', 'r44266240', 'r44266496', 'r44266752', 'r44267008', 'r44267264', 'r44267520', 'r44267776', 'r44268032', 'r44268288', 'r44268544', 'r44268800', 'r44269056', 'r44269312', 'r44269568', 'r44269824', 'r44270080', , 'r45493248', 'r45493504', 'r45493760', 'r45494016', 'r45494272', 'r45494528', 'r45494784', 'r45495040', 'r45495296', 'r45495552',
;aorname_ch2 = [ 'r42532864',  'r42533120']
  basedir = '/Users/jkrick/irac_warm/calstars/' 
  utmjd_center =   56701.677D ; junk
  transit_duration=69.6  ; junk
  period = 0.813475       ;junk
  intended_phase = 0.5
  exptime = 0.1
  stareaor = 1E5
  plot_norm =0.3852
  plot_corrnorm = 0.3852

   values=list(ra_ref, dec_ref, aorname_ch2, aorname_ch1, basedir, chname, utmjd_center, transit_duration, period, intended_phase, exptime, mask, maskreg, stareaor, plot_norm, plot_corrnorm)
  planetinfo[planetname] = HASH(keys, values)



;  print, planetinfo.keys()
;  print, planetinfo['wasp13'].keys()
;  print, 'plan', planetinfo['wasp15', 'aorname']
  
  return, planetinfo

end
