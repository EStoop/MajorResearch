function find_chname, dirname
  ;;assume neither are there to start
  ch2 = 0
  ch1 = 0

  ch2 = file_test(strcompress(dirname + '/ch2',/remove_all), /Directory)
  ch1 = file_test(strcompress(dirname + '/ch1',/remove_all), /Directory)

  if ch1 eq 1 then chname = 'ch1'
  if ch2 eq 1 then chname = 'ch2'
  if ch1 eq 1 and ch2 eq 1 then chname = ['ch1', 'ch2']
  
  return, chname
end


;  cd, dirname
;  lscommand = 'ls > dirlist.txt'
;  spawn, lscommand;;

;  readcol, 'dirlist.txt', filelist, format = '(A)'
;  for i = 0, n_elements(filelist) - 1 do begin
;     print, 'filelist', filelist
;     if filelist(i) eq 'ch2' then ch2 = 1
;     if filelist(i) eq 'ch1' then ch1 = 1
;  endfor
