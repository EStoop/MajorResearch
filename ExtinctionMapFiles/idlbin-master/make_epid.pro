pro make_epid

;this code takes a skydarks_XX file and makes an EPID file
;for use in processing the dark files

base_dir = '/Users/jkrick/irac_darks/pipeline/'
campaign_name = 'pc83'

;read in the first line to get the directory name and pertinent info
readcol, strcompress(base_dir +'skydarks_' + campaign_name,/remove_all), dir_name , camp_name, junk, junk2,  process_num, format = 'A, A, A,A,A'
process = strcompress('S'+process_num,/remove_all)

;then read in the rest of the lines with the info
readcol, strcompress(base_dir +'skydarks_' + campaign_name,/remove_all), aor, good, name,  format = 'L, A, A'


;open the output EPID file
openw, outlun,strcompress(base_dir +'EPID_' + campaign_name + '_' + process,/remove_all), /get_lun
printf, outlun, '     Id            Type    Instrument        Channel       Status'

;for each aor, 
for i = 0, n_elements(aor) - 1 do begin
   cd,strcompress( '~/iracdata/flight/IWIC/'+ dir_name +'/cal/00'+string(aor(i)),/remove_all)
   command = 'ls *.fits > '+strcompress(base_dir+'fits_list.txt',/remove_all)
   spawn, command

   ;read in the fits file headers
   readcol, strcompress(base_dir + 'fits_list.txt',/remove_all), fitsname, format = 'A'

   ;print to the EPID file the appropriate EPID's
    ;take into account that some of the
    ;darks will be bad as indicated by the
    ;second column of the skydarks file
   for n = 0, n_elements(fitsname) - 1 do begin
      header = headfits(fitsname(n)) 
      ch = fxpar(header, 'CHNLNUM')

      print, 'good, ch epid', good(i), ch,  fxpar(header, 'EPID')
      if good(i) eq 'all' then  printf, outlun, fxpar(header, 'EPID'), '    m           IRAC      ', ch , '             2048'
      if good(i) eq '1' and ch eq 1 then printf, outlun, fxpar(header, 'EPID'), '    m           IRAC      ', ch , '             2048'
      if good(i) eq '2' and ch eq 2 then printf, outlun, fxpar(header, 'EPID'), '    m           IRAC      ', ch , '             2048'
     
   endfor
   printf, outlun, "  "
endfor
close, outlun
free_lun, outlun
   
end
