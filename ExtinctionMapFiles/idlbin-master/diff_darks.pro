pro diff_darks
  
;this code makes a jpeg that is the difference between 2 sets of darks taken during one warm campaign.  This is done for all exptimes.

base_dir = '/Users/jkrick/irac_darks/pipeline/'
campaign_name = 'pc83'

;read in the first line to get the directory name and pertinent info
readcol, strcompress(base_dir +'skydarks_' + campaign_name,/remove_all), dir_name , camp_name, junk, junk2,  process_num, format = 'A, A, A,A,A',/silent
process = strcompress('S'+process_num,/remove_all)

;then read in the rest of the lines with the info
readcol, strcompress(base_dir +'skydarks_' + campaign_name,/remove_all), aor, good, name,  format = 'L, A, A',/silent

;print, 'working on ', dir_name

;for each matching set of aors, not including the subarrays
for i = 4,  13 do begin
;which other directory has the matching exptime?
   j = i + 14

   ;print, 'working on ', aor(i), aor(j)

;open up first set of images (a) for this exptime.
   cd,strcompress( '~/iracdata/flight/IWIC/'+ dir_name +'/cal/00'+string(aor(i)),/remove_all)
   command = 'ls *.fits > '+strcompress(base_dir+'fits_a_list.txt',/remove_all)
   spawn, command

;open up second set of images (b) for this exptime.
   cd,strcompress( '~/iracdata/flight/IWIC/'+ dir_name +'/cal/00'+string(aor(j)),/remove_all)
   command = 'ls *.fits > '+strcompress(base_dir+'fits_b_list.txt',/remove_all)
   spawn, command

   ;read in the fits file names
   readcol, strcompress(base_dir + 'fits_a_list.txt',/remove_all), fitsname_a, format = 'A',/silent
   readcol, strcompress(base_dir + 'fits_b_list.txt',/remove_all), fitsname_b, format = 'A',/silent

   ;for each file in that aor, 2 if full array, non-hdr, more for hdr.
    for n = 0, n_elements(fitsname_a) - 1 do begin

       filenm_a = strcompress( '~/iracdata/flight/IWIC/'+ dir_name +'/cal/00'+string(aor(i)) + '/'+ fitsname_a(n),/remove_all)
      ; header_a = headfits(filenm_a) 
       fits_read, filenm_a, data_a, header_a
       e_a = fxpar(header_a, 'FRAMTIME')
       ch_a = fxpar(header_a, 'CHNLNUM')

       filenm_b = strcompress( '~/iracdata/flight/IWIC/'+ dir_name +'/cal/00'+string(aor(j)) + '/'+ fitsname_b(n),/remove_all)
       ;header_b = headfits(filenm_b) 
       fits_read, filenm_b, data_b, header_b
       e_b = fxpar(header_b, 'FRAMTIME')
       ch_b = fxpar(header_b, 'CHNLNUM')

       ;print, 'differenconfg', e_a, e_b, ch_a, ch_b
       
       diff = data_a - data_b
       ;print,   strcompress(base_dir+campaign_name +'_diff_' + string(ch_a) +'_'+ string(e_a) + '.fits',/remove_all)
       fits_write, strcompress(base_dir+campaign_name +'_diff_' + string(ch_a) +'_'+ string(e_a) + '.fits',/remove_all), diff, header_a
       
   endfor
endfor  
 

end
