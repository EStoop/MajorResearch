pro diff_headers, aor, channel

dir = '~/iracdata/FINAL_CHECK/IRAC012200/'
ndir = dir + 'S18.18v2/'
odir = dir + 'S18.7/'
saor = strn(aor)
sch = strn(channel)
suffix = 'r' + saor + '/ch' + sch + '/bcd/'
sstr = ndir + suffix
xstr = odir + suffix

files = file_search(sstr + '/SPITZER_I' + sch + '_' + saor + '_*_cbcd.fits', COUNT=nfiles)
xfiles = file_search(xstr + '/SPITZER_I' + sch + '_' + saor + '_*_cbcd.fits', COUNT=xcount)
if (xcount eq nfiles) then begin
   for j = 10, 10 do begin; nfiles-1 do begin
           
      nbcdhead = headfits(files[j]) ;
      obcdhead = headfits(xfiles[j])

      nh = sstr + 'nbcdhead.txt'
      openw, noutlun, nh, /get_lun
      printf, noutlun, nbcdhead
      close, noutlun
      free_lun, noutlun

      oh = sstr + 'obcdhead.txt'
      openw, ooutlun, oh, /get_lun
      printf, ooutlun, obcdhead
      close, ooutlun
      free_lun, ooutlun
    
      command = 'diff -l ' + nh + ' ' + oh
      print, 'running diff'
      spawn, command
   endfor
endif




end
