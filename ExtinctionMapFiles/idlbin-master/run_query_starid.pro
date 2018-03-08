pro run_query_starid

  cd, '/Users/jkrick/external/irac_warm/trending/r63202816/ch2/bcd/'

  readcol,'bcdlist.txt',fitsname, format = 'A', /silent
  
  
  header = headfits(fitsname(0))
  ra_rqst = sxpar(header, 'RA_RQST') ;need to use these so that I know when there is a blank channel
  dec_rqst = sxpar(header, 'DEC_RQST')
  naxishere = sxpar(header, 'NAXIS')
  date_obs = sxpar(header, 'DATE_OBS')
  bmjd0 = sxpar(header, 'BMJD_OBS')
  year_obs = float(date_obs.substring(0,3))
  ra = ra_rqst
  dec = dec_rqst
  starname = Query_starid_v2( header, ra = ra, dec = dec ,/Verbose)
  print, 'starname: ', starname, ra, dec
  
  ;;result = get_exoplanet_info(ra, dec)
  ;;print, result
  ;;how many objects were found?
 end
