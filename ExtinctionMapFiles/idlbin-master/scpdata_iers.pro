function scpdata_iers, aorname, campaign_name
  ;;copy data from sscdev100 archive to my machine
  
  ;;start with making the directory to house it
  basedir = '/Volumes/external/'
  
  localname =strcompress(basedir +'irac_warm/trending/r' + string(aorname),/remove_all)
;;  print, 'localname ', localname
  spawn, 'mkdir '+ localname
  cd, localname
  if file_test('*') lt 1 then begin
     print, 'geting data'
     remotename = strcompress('jkrick@sscdev100.ipac.caltech.edu:/sha/archive/proc/'+campaign_name + '/r' + string(aorname) ,/remove_all)
     scpcommand = 'scp -qr '+ remotename + '  .'
     print, 'scp command', scpcommand
     spawn, scpcommand
  endif

  ;;now to figure out which channels there are
  spawn, "ls > '../chname.txt'"
  readcol, '../chname.txt', junkstring, format = '(A)'
  fullstr = ''
  for i = 0, n_elements(junkstring) - 1 do fullstr = fullstr + junkstring(i)
  if fullstr.contains('ch1') then begin
     chname = 'ch1'
     cd, chname
     spawn, 'rm -rf cal/'
     spawn, 'rm -rf pbcd/'
     cd, 'bcd'
     spawn, "rm -rf *cbcd.fits"
     spawn, "rm -rf *cbunc.fits"
     spawn, "rm -rf *msk.fits"

  endif
  ;; preference to reduce ch2 if both exist
  if fullstr.contains('ch2') then begin
     chname = 'ch2'             
     cd, chname
     spawn, 'rm -rf cal/'
     spawn, 'rm -rf pbcd/'
     cd, 'bcd'
     spawn, "rm -rf *cbcd.fits"
     spawn, "rm -rf *cbunc.fits"
     spawn, "rm -rf *msk.fits"
 endif
  

  
   
  if file_test('*dce.fits') lt 1 then begin
     print, 'getting dce'
     remoteraw = strcompress('jkrick@sscdev100.ipac.caltech.edu:/stage/sha/archive/raw/campaign/'+campaign_name + '/r' + string(aorname) + '/' + chname + '/raw/\*_dce.fits',/remove_all)
     scpraw = 'scp -q '+ remoteraw + '  .'
     print, 'scp dce command', scpraw
     spawn, scpraw
  endif
  
  return, 0
end
