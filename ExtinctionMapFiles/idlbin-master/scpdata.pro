function scpdata, aorname, campaign_name, chname
  ;;copy data from sscdev100 archive to my machine
  
  ;;start with making the directory to house it
  basedir = '/Volumes/external/'
  
  localname =strcompress(basedir +'irac_warm/trending/r' + string(aorname),/remove_all)
;;  print, 'localname ', localname
  spawn, 'mkdir '+ localname
  cd, localname
  spawn, 'mkdir '+ chname
  cd, chname
  spawn, 'mkdir bcd'
  cd, 'bcd'
  if file_test('*bcd.fits') lt 1 then begin
     print, 'geting bcd'
     remotename = strcompress('jkrick@sscdev100.ipac.caltech.edu:/sha/archive/proc/'+campaign_name + '/r' + string(aorname) + '/' + chname + '/bcd/\*_bcd.fits',/remove_all)
     scpcommand = 'scp -q '+ remotename + '  .'
     print, 'scp command', scpcommand
     spawn, scpcommand
  endif
  
  if file_test('*bunc.fits') lt 1 then begin
     print, 'getting bunc'
     remotebunc = strcompress('jkrick@sscdev100.ipac.caltech.edu:/sha/archive/proc/'+campaign_name + '/r' + string(aorname) + '/' + chname + '/bcd/\*_bunc.fits',/remove_all)
     scpbunc = 'scp -q '+ remotebunc + '  .'
     spawn, scpbunc
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
