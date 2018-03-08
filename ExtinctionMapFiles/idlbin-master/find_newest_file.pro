function find_newest_file
  dirname = '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/'
  cd, dirname
  command  ="find . -name 'track*.sav' > ../filelist.txt"
  spawn, command
  
  readcol, '../filelist.txt',savename, format = 'A', /silent
  print, 'nelements savename', n_elements(savename)
  if n_elements(savename) gt 0 then begin
     mtime = file_modtime(strcompress(dirname + savename))
     junk = max(mtime, m)
     newest_save = dirname + savename(m).substring(2)
     print, 'newest: ',newest_save
     return, newest_save
  endif else begin
     return, 'thereisnofile.sav'
  endelse
  
end
