function quick_outlier, array, nsig
  m = median(array)
  st = stddev(array, /nan)
  print, 'criteria', m, nsig*st
  good = where(array lt (m + nsig*st) and array gt (m - nsig*st), complement = bad)
  return, bad
end



pro save_track_centroids

  ;;this code writes out a csv file to use in python to attempt pretty plots
   restore,  '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/plot_track_centroids_test.sav'
   laorname =  long(aorname)
   xjd = double(xjd)
   print, 'laorname', laorname[0:10]
   ;;get rid of nans, these are likely the short AORs
   good = where(finite(npmean) gt 0) ;bunch of nans in this array
   npmean = npmean(good)
   laorname=laorname(good)
   exptimearr=exptimearr(good)
   xjd=xjd(good)
   startyear=startyear(good)
   startmonth = startmonth(good)
   pa=pa(good)
   dpa=dpa(good)
   sigmax=sigmax(good)
   sigmay=sigmay(good)
   xdrift=xdrift(good)
   ydrift=ydrift(good)
   short_drift=short_drift(good)
   slope_drift=slope_drift(good)
   pkperiod=pkperiod(good)
   pkstrength=pkstrength(good)
   npunc=npunc(good)   
   avgxcen = avgxcen(good)
   avgycen = avgycen(good)
   obsdur = obsdur(good)
   pitch_earth = pitch_earth(good)
   rms_slope_pld=rms_slope_pld(good)
   rms_absdev_pld= rms_absdev_pld(good)
   rms_slope_raw=rms_slope_raw(good)
   rms_absdev_raw= rms_absdev_raw(good)
   chname = chname(good)
   
   ;;remove some outliers
   b = where(xdrift gt 0.3)  ;;can't figure out otherwise how to narrow the range, there are only a few out there
   xdrift(b) = alog10(-1)
   b = where(ydrift gt 0.3)
   ydrift(b) = alog10(-1)
   b = where(pkstrength lt 0)    ;;shouldn't be negative values in these arrays
   pkstrength(b) = alog10(-1)
   b = where(sigmay gt 1)
   sigmay(b) = alog10(-1)
   b = where(sigmax gt 1)
   sigmax(b) = alog10(-1)

   
   ;;three sigma on npmean, ydrift, xdrift, sigmax, sigmay
   ;;write my own quick outlier rejection since I can;t get any of the
   ;;IDL ones to work.

   bad = quick_outlier(npmean, 3.0)
   npmean(bad) = alog10(-1)
   print, 'n_elements npmean outlier', n_elements(bad)
   
   bad = quick_outlier(sigmax, 3.0)
   sigmax(bad) = alog10(-1)
   print, 'n_elements sigmax outlier', n_elements(bad)

   bad = quick_outlier(sigmay, 3.0)
   sigmay(bad) = alog10(-1)
   print, 'n_elements sigmay outlier', n_elements(bad)

   bad = quick_outlier(xdrift, 3.0)
   xdrift(bad) = alog10(-1)
   print, 'n_elements xdrift outlier', n_elements(bad)

   bad = quick_outlier(ydrift, 3.0)
   ydrift(bad) = alog10(-1)
   print, 'n_elements ydrift outlier', n_elements(bad)



   ;;force similar exptimes into the same bin
   a = where(exptimearr eq 23.6)
   exptimearr(a) = 26.8
   a = where(exptimearr eq 93.6)
   exptimearr(a) = 96.8

   print, 'laorname before array concatenation', [laorname[0:10]]
   help, laorname

   ;;look for certain aors
   g = where(xjd gt  2457828.5, ng)
   print, 'ng', ng
   print, xjd(g)

   ;;want to remove either x or ydrift = -0.05
   baddrift = where(xdrift lt-0.049995 and xdrift gt -0.050001, nbaddrift)
   xdrift(baddrift) = alog10(-1)
   baddrift = where(ydrift lt-0.049995 and ydrift gt -0.050001, nbaddrift)
   ydrift(baddrift) = alog10(-1)
   
   ;; chname can't be a string to use array concatenation
   ;;change to a number instead of figuring out diferent data
   ;;structures now that I know this one works
   chname_dble = double(chname.Substring(2))
   
   ;;array concatenation fun
   data =[[laorname],[exptimearr],[xjd], [startyear],[startmonth], [pa],[dpa], [sigmax], [sigmay], [xdrift], [ydrift],[short_drift], [slope_drift], [pkperiod], [pkstrength], [npmean], [npunc] , [avgxcen], [avgycen],[obsdur*24.],[pitch_earth],[rms_slope_raw],[rms_absdev_raw],[rms_slope_pld], [rms_absdev_pld],[chname_dble]]

   print, 'testong', data[0]
   help, data, exptimearr, xdrift, pkstrength, laorname, chname
   datac = transpose(data)

   ;;write out the csv file for python plotting
   h = ['AORNAME','EXPTIME','JD', 'STARTYEAR','STARTMONTH','PITCH_ANGLE', 'DELTA_PA', 'SIGMAX', 'SIGMAY', 'XDRIFT', 'YDRIFT', 'SHORT_DRIFT', 'SLOPE_DRIFT', 'PKPERIOD', 'PKSTRENGTH', 'NPMEAN', 'NPUNC','AVGXCEN','AVGYCEN','OBSDUR','PITCH_EARTH','RMS_SLOPE_RAW', 'RMS_ABSDEV_RAW','RMS_SLOPE_PLD','RMS_ABSDEV_PLD','CHNAME']
   write_csv, '~/pybin/track_centroids.csv', datac, $
              header = h


end


