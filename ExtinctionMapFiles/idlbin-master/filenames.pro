pro filenames, dir_name;, aor_list

;dir_name = '/Users/jkrick/iwic/IWICA/IRAC016300'

close, /all

exptime = [0.1,0.4,2,6,12,30, 100,200]
cd, dir_name
command = 'ls raw | grep -v raw > aorlist.txt'
;spawn, command
readcol, 'aorlist.txt', aorname, format="A",/silent

;open file lists for all three positions together

openw,90, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_200s.txt')
openw,89, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_100s.txt')
openw,88, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_200_64s.txt')
openw,87, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_100_32s.txt')
openw,86, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_30s.txt')
openw,85, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_12s.txt')
openw,84, strcompress(dir_name+ '/raw/' +  'rawlist_ch1_6s.txt')
openw,83, strcompress(dir_name+ '/raw/' +  'rawlist_ch1_2s.txt')
openw,91, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_0pt6s.txt')
openw,82, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_0pt4s.txt')
openw,81, strcompress(dir_name+ '/raw/' + 'rawlist_ch1_0pt1s.txt')
openw,80, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_200s.txt')
openw,79, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_100s.txt')
openw,78, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_200_64s.txt')
openw,77, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_100_32s.txt')
openw,76, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_30s.txt')
openw,75, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_12s.txt')
openw,74, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_6s.txt')
openw,73, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_2s.txt')
openw,72, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_0pt4s.txt')
openw,92, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_0pt6s.txt')
openw,71, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch1_0pt1s.txt')
openw,70, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_200s.txt')
openw,69, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_100s.txt')
openw,68, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_200_64s.txt')
openw,67, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_100_32s.txt')
openw,66, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_30s.txt')
openw,65, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_12s.txt')
openw,64, strcompress(dir_name+ '/raw/' +  'rawlist_ch2_6s.txt')
openw,63, strcompress(dir_name+ '/raw/' +  'rawlist_ch2_2s.txt')
openw,93, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_0pt6s.txt')
openw,62, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_0pt4s.txt')
openw,61, strcompress(dir_name+ '/raw/' + 'rawlist_ch2_0pt1s.txt')
openw,60, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_200s.txt')
openw,59, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_100s.txt')
openw,58, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_200_64s.txt')
openw,57, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_100_32s.txt')
openw,56, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_30s.txt')
openw,55, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_12s.txt')
openw,54, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_6s.txt')
openw,53, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_2s.txt')
openw,94, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_0pt6s.txt')
openw,52, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_0pt4s.txt')
openw,51, strcompress(dir_name+ '/bcd/' + 'bcdlist_ch2_0pt1s.txt')



;for each of the AORs (3)
for n = 0,  n_elements(aorname) - 1 do begin


   print, 'working on aor', aorname(n)
   ;change to that AOR directory
;make lists by channel number for raw and bcd
   CD, strcompress(dir_name + '/raw/' + aorname(n))

   command1 =  ' find . -name "*.fits" >  raw_list.txt'
   command2 = 'cat raw_list.txt | grep IRAC.1. > raw_ch1.txt'
   command3 = 'cat raw_list.txt | grep IRAC.2. > raw_ch2.txt'

   commands = [command1, command2, command3]
   for i = 0, n_elements(commands) -1 do spawn, commands(i)

;now bcd
   CD, strcompress(dir_name + '/bcd/' + aorname(n))
   command1 =  ' find . -name "*bcd_fp.fits" >  bcd_list.txt'
   command2 = 'cat bcd_list.txt | grep IRAC.1. > bcd_ch1.txt'
   command3 = 'cat bcd_list.txt | grep IRAC.2. > bcd_ch2.txt'

   commands = [command1, command2, command3]
   for i = 0, n_elements(commands) -1 do spawn, commands(i)

;now need to make lists by exposure time
;which means reading in the headers of these fits files.  all of them
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;ch1
;start with the raws
   readcol, strcompress(dir_name + '/raw/'+ aorname(n) + '/raw_ch1.txt'), rawch1list, format="A", /silent
   count200 = 0
   list200 = strarr(200)
   count100 = 0
   list100 = strarr(200)
   count200_64 = 0
   list200_64 = strarr(200)
   count100_32 = 0
   list100_32 = strarr(200)
   count30 = 0
   list30 = strarr(200)
   count12 = 0
   list12 = strarr(200)
   count6 = 0
   list6 = strarr(200)
   count2 = 0
   list2 = strarr(200)
   count0pt6 = 0
   list0pt6 = strarr(200)
   count0pt4 = 0
   list0pt4 = strarr(200)
   count0pt1 = 0
   list0pt1 = strarr(200)

   CD, strcompress(dir_name + '/raw/' + aorname(n))

   for j = 0, n_elements(rawch1list) - 1 do begin
      rawheader = headfits(rawch1list(j)) ;
      fowlnum = fxpar(rawheader, 'A0614D00')
      waitper = fxpar(rawheader, 'A0615D00')

      if fowlnum eq 32 and waitper eq 936 then begin
         list200(count200) = rawch1list(j)
         count200 = count200+1
      endif
      if fowlnum eq 16 and waitper eq 468 then begin
         list100(count100) = rawch1list(j)
         count100 = count100+1
      endif
      if fowlnum eq 64 and waitper eq 872 then begin
         list200_64(count200_64) = rawch1list(j)
         count200_64 = count200_64+1
      endif
      if fowlnum eq 32 and waitper eq 436 then begin
         list100_32(count100_32) = rawch1list(j)
         count100_32 = count100_32+1
      endif
      if fowlnum eq 16 and waitper eq 118 then begin
         list30(count30) = rawch1list(j)
         count30 = count30+1
      endif
     if fowlnum eq 8 and waitper eq 44 then begin
         list12(count12) = rawch1list(j)
         count12 = count12+1
      endif
      if fowlnum eq 8 and waitper ne 44 then begin
         list6(count6) = rawch1list(j)
         count6 = count6+1
      endif
      if fowlnum eq 4 and waitper eq 2 then begin
         list2(count2) = rawch1list(j)
         count2 = count2+1
      endif
      if fowlnum eq 1 and waitper eq 1 then begin
         list0pt6(count0pt6) = rawch1list(j)
         count0pt6 = count0pt6+1
      endif
      if fowlnum eq 1 and waitper eq 0 then begin
         list0pt4(count0pt4) = rawch1list(j)
         count0pt4 = count0pt4+1
      endif
      if fowlnum eq 2 and waitper eq 6 then begin
         list0pt1(count0pt1) = rawch1list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor


   openw, outlun200, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,90, strcompress(aorname(n) + strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 89, strcompress(aorname(n) + strmid(list100(i),1))
   endfor

  openw, outlun200_64, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_200_64s.txt'),/GET_LUN
   for i = 0, count200_64 - 1 do begin
      printf, outlun200_64, list200_64(i)
      printf,88, strcompress(aorname(n) + strmid(list200_64(i),1))
   endfor

   openw, outlun100_32, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_100_32s.txt'),/GET_LUN
   for i = 0, count100_32 - 1 do begin
      printf, outlun100_32, list100_32(i)
      printf, 87, strcompress(aorname(n) + strmid(list100_32(i),1))
   endfor

   openw, outlun30, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_30s.txt'),/GET_LUN
   for i = 0, count30 - 1 do begin
      printf, outlun30, list30(i)
      printf, 86, strcompress(aorname(n) + strmid(list30(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 85, strcompress(aorname(n) + strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 84, strcompress(aorname(n) + strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 83, strcompress(aorname(n) + strmid(list2(i),1))
   endfor

   openw, outlun0pt6, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_0pt6s.txt'),/GET_LUN
   for i = 0, count0pt6 - 1 do begin
      printf, outlun0pt6, list0pt6(i)
      printf, 91, strcompress(aorname(n) + strmid(list0pt6(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 82, strcompress(aorname(n) + strmid(list0pt4(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch1_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 81, strcompress(aorname(n) + strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
   close, outlun200_64
   close, outlun100_32
   close, outlun30
   close, outlun6
   close, outlun2
   close, outlun0pt4
   close, outlun0pt6
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun200_64
   free_lun, outlun100_32
   free_lun, outlun30
   free_lun, outlun12
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt4
   free_lun, outlun0pt6
   free_lun, outlun0pt1
   

;----------------------------------------------------------------------------------
;----------------------------------------------------------------------------------
;could fake the bcd's
;faster just to copy everything from above

;go on the bcd's
   print, strcompress(dir_name + '/bcd/'+ aorname(n) + '/bcd_ch1.txt')
   readcol, strcompress(dir_name + '/bcd/'+ aorname(n) + '/bcd_ch1.txt'), bcdch1list, format="A", /silent
   count200 = 0
   list200 = strarr(200)
   count100 = 0
   list100 = strarr(200)
   count12 = 0
   list12 = strarr(200)
   count200_64 = 0
   list200_64 = strarr(200)
   count100_32 = 0
   list100_32 = strarr(200)
   count30 = 0
   list30 = strarr(200)
   count6 = 0
   list6 = strarr(200)
   count2 = 0
   list2 = strarr(200)
   count0pt6 = 0
   list0pt6 = strarr(200)
   count0pt4 = 0
   list0pt4 = strarr(200)
   count0pt1 = 0
   list0pt1 = strarr(200)

   CD, strcompress(dir_name + '/bcd/' + aorname(n))

   for j = 0, n_elements(bcdch1list) - 1 do begin
      bcdheader = headfits(bcdch1list(j)) ;
      frametime = fxpar(bcdheader, 'framtime')
      fowlnum = fxpar(bcdheader, 'AFOWLNUM')
      waitper = fxpar(bcdheader, 'AWAITPER')

      if frametime eq 200 and fowlnum eq 32 then begin
         list200(count200) = bcdch1list(j)
         count200 = count200+1
      endif
      if frametime eq 100 and fowlnum eq 16 then begin
         list100(count100) = bcdch1list(j)
         count100 = count100+1
      endif
      if fowlnum eq 64 and waitper eq 872 then begin
         list200_64(count200_64) = bcdch1list(j)
         count200_64 = count200_64+1
      endif
      if fowlnum eq 32 and waitper eq 436 then begin
         list100_32(count100_32) = bcdch1list(j)
         count100_32 = count100_32+1
      endif
      if fowlnum eq 16 and waitper eq 118 then begin
         list30(count30) = bcdch1list(j)
         count30 = count30+1
      endif
      if frametime eq 12 then begin
         list12(count12) = bcdch1list(j)
         count12 = count12+1
      endif
      if frametime eq 6 then begin
         list6(count6) = bcdch1list(j)
         count6 = count6+1
      endif
      if frametime eq 2 then begin
         list2(count2) = bcdch1list(j)
         count2 = count2+1
      endif
      if frametime eq 0.6 then begin
         list0pt6(count0pt6) = bcdch1list(j)
         count0pt6 = count0pt6+1
      endif
      if frametime eq 0.4 then begin
         list0pt4(count0pt4) = bcdch1list(j)
         count0pt4 = count0pt4+1
      endif
      if frametime eq 0.1 then begin
         list0pt1(count0pt1) = bcdch1list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor

  
 openw, outlun200, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,80, strcompress(aorname(n) + strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_100s.txt'),/GET_LUN
   print, 'aorname(n), count100', aorname(n), count100
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 79, strcompress(aorname(n) + strmid(list100(i),1))
   endfor
 openw, outlun200_64, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_200_64s.txt'),/GET_LUN
   for i = 0, count200_64 - 1 do begin
      printf, outlun200_64, list200_64(i)
      printf,78, strcompress(aorname(n) + strmid(list200_64(i),1))
   endfor

   openw, outlun100_32, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_100_32s.txt'),/GET_LUN
   for i = 0, count100_32 - 1 do begin
      printf, outlun100_32, list100_32(i)
      printf, 77, strcompress(aorname(n) + strmid(list100_32(i),1))
   endfor

   openw, outlun30, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_30s.txt'),/GET_LUN
   for i = 0, count30 - 1 do begin
      printf, outlun30, list30(i)
      printf, 76, strcompress(aorname(n) + strmid(list30(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 75, strcompress(aorname(n) + strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 74, strcompress(aorname(n) + strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 73, strcompress(aorname(n) + strmid(list2(i),1))
   endfor

   openw, outlun0pt6, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_0pt6s.txt'),/GET_LUN
   for i = 0, count0pt6 - 1 do begin
      printf, outlun0pt6, list0pt6(i)
      printf, 92, strcompress(aorname(n) + strmid(list0pt6(i),1))
   endfor
      
   openw, outlun0pt4, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 72, strcompress(aorname(n) + strmid(list0pt4(i),1))
   endfor

   openw, outlun0pt1, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch1_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 71, strcompress(aorname(n) + strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
 close, outlun200_64
   close, outlun100_32
   close, outlun30
   close, outlun6
   close, outlun2
   close, outlun0pt6
   close, outlun0pt4
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun12
   free_lun, outlun200_64
   free_lun, outlun100_32
   free_lun, outlun30
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt6
   free_lun, outlun0pt4
   free_lun, outlun0pt1


;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;ch2
;start with the raws
   readcol, strcompress(dir_name + '/raw/'+ aorname(n) + '/raw_ch2.txt'), rawch2list, format="A", /silent
  count200 = 0
   list200 = strarr(200)
   count100 = 0
   list100 = strarr(200)
   count200_64 = 0
   list200_64 = strarr(200)
   count100_32 = 0
   list100_32 = strarr(200)
   count30 = 0
   list30 = strarr(200)
   count12 = 0
   list12 = strarr(200)
   count6 = 0
   list6 = strarr(200)
   count2 = 0
   list2 = strarr(200)
   count0pt6 = 0
   list0pt6 = strarr(200)
   count0pt4 = 0
   list0pt4 = strarr(200)
   count0pt1 = 0
   list0pt1 = strarr(200)

   CD, strcompress(dir_name + '/raw/' + aorname(n))

   for j = 0, n_elements(rawch2list) - 1 do begin
      rawheader = headfits(rawch2list(j)) ;
      fowlnum = fxpar(rawheader, 'A0614D00')
      waitper = fxpar(rawheader, 'A0615D00')

      if fowlnum eq 32 and waitper eq 936 then begin
         list200(count200) = rawch2list(j)
         count200 = count200+1
      endif
      if fowlnum eq 16 and waitper eq 468 then begin
         list100(count100) = rawch2list(j)
         count100 = count100+1
      endif
      if fowlnum eq 64 and waitper eq 872 then begin
         list200_64(count200_64) = rawch2list(j)
         count200_64 = count200_64+1
      endif
      if fowlnum eq 32 and waitper eq 436 then begin
         list100_32(count100_32) = rawch2list(j)
         count100_32 = count100_32+1
      endif
      if fowlnum eq 16 and waitper eq 118 then begin
         list30(count30) = rawch2list(j)
         count30 = count30+1
      endif
     if fowlnum eq 8 and waitper eq 44 then begin
         list12(count12) = rawch2list(j)
         count12 = count12+1
      endif
      if fowlnum eq 8 and waitper ne 44 then begin
         list6(count6) = rawch2list(j)
         count6 = count6+1
      endif
      if fowlnum eq 4 and waitper eq 2 then begin
         list2(count2) = rawch2list(j)
         count2 = count2+1
      endif
      if fowlnum eq 1 and waitper eq 1 then begin
         list0pt6(count0pt6) = rawch2list(j)
         count0pt6 = count0pt6+1
      endif
      if fowlnum eq 1 and waitper eq 0 then begin
         list0pt4(count0pt4) = rawch2list(j)
         count0pt4 = count0pt4+1
      endif
      if fowlnum eq 2 and waitper eq 6 then begin
         list0pt1(count0pt1) = rawch2list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor


   openw, outlun200, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,70, strcompress(aorname(n) + strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 69, strcompress(aorname(n) + strmid(list100(i),1))
   endfor

  openw, outlun200_64, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_200_64s.txt'),/GET_LUN
   for i = 0, count200_64 - 1 do begin
      printf, outlun200_64, list200_64(i)
      printf,68, strcompress(aorname(n) + strmid(list200_64(i),1))
   endfor

   openw, outlun100_32, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_100_32s.txt'),/GET_LUN
   for i = 0, count100_32 - 1 do begin
      printf, outlun100_32, list100_32(i)
      printf, 67, strcompress(aorname(n) + strmid(list100_32(i),1))
   endfor

   openw, outlun30, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_30s.txt'),/GET_LUN
   for i = 0, count30 - 1 do begin
      printf, outlun30, list30(i)
      printf, 66, strcompress(aorname(n) + strmid(list30(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 65, strcompress(aorname(n) + strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 64, strcompress(aorname(n) + strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 63, strcompress(aorname(n) + strmid(list2(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 62, strcompress(aorname(n) + strmid(list0pt4(i),1))
   endfor

   openw, outlun0pt6, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_0pt6s.txt'),/GET_LUN
   for i = 0, count0pt6 - 1 do begin
      printf, outlun0pt6, list0pt6(i)
      printf, 93, strcompress(aorname(n) + strmid(list0pt6(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/raw/' + aorname(n)+ '/rawlist_ch2_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 61, strcompress(aorname(n) + strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
   close, outlun200_64
   close, outlun100_32
   close, outlun30
   close, outlun6
   close, outlun2
   close, outlun0pt4
   close, outlun0pt6
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun200_64
   free_lun, outlun100_32
   free_lun, outlun30
   free_lun, outlun12
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt6
   free_lun, outlun0pt4
   free_lun, outlun0pt1
   


;----------------------------------------------------------------------------------
;----------------------------------------------------------------------------------
;could fake the bcd's
;faster just to copy everything from above

;go on the bcd's
   readcol, strcompress(dir_name + '/bcd/'+ aorname(n) + '/bcd_ch2.txt'), bcdch2list, format="A", /silent
   count200 = 0
   list200 = strarr(200)
   count100 = 0
   list100 = strarr(200)
   count12 = 0
   list12 = strarr(200)
   count200_64 = 0
   list200_64 = strarr(200)
   count100_32 = 0
   list100_32 = strarr(200)
   count30 = 0
   list30 = strarr(200)
   count6 = 0
   list6 = strarr(200)
   count2 = 0
   list2 = strarr(200)
   count0pt6 = 0
   list0pt6 = strarr(200)
  count0pt4 = 0
   list0pt4 = strarr(200)
   count0pt1 = 0
   list0pt1 = strarr(200)

   CD, strcompress(dir_name + '/bcd/' + aorname(n))

   for j = 0, n_elements(bcdch2list) - 1 do begin
      bcdheader = headfits(bcdch2list(j)) ;
      frametime = fxpar(bcdheader, 'framtime')
      fowlnum = fxpar(bcdheader, 'AFOWLNUM')
      waitper = fxpar(bcdheader, 'AWAITPER')

      if frametime eq 200 and fowlnum eq 32 then begin
         list200(count200) = bcdch2list(j)
         count200 = count200+1
      endif
      if frametime eq 100 and fowlnum eq 16 then begin
         list100(count100) = bcdch2list(j)
         count100 = count100+1
      endif
      if fowlnum eq 64 and waitper eq 872 then begin
         list200_64(count200_64) = bcdch2list(j)
         count200_64 = count200_64+1
      endif
      if fowlnum eq 32 and waitper eq 436 then begin
         list100_32(count100_32) = bcdch2list(j)
         count100_32 = count100_32+1
      endif
      if fowlnum eq 16 and waitper eq 118 then begin
         list30(count30) = bcdch2list(j)
         count30 = count30+1
      endif
      if frametime eq 12 then begin
         list12(count12) = bcdch2list(j)
         count12 = count12+1
      endif
      if frametime eq 6 then begin
         list6(count6) = bcdch2list(j)
         count6 = count6+1
      endif
      if frametime eq 2 then begin
         list2(count2) = bcdch2list(j)
         count2 = count2+1
      endif
      if frametime eq 0.6 then begin
         list0pt6(count0pt6) = bcdch2list(j)
         count0pt6 = count0pt6+1
      endif
      if frametime eq 0.4 then begin
         list0pt4(count0pt4) = bcdch2list(j)
         count0pt4 = count0pt4+1
      endif
      if frametime eq 0.1 then begin
         list0pt1(count0pt1) = bcdch2list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor

  
 openw, outlun200, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,60, strcompress(aorname(n) + strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 59, strcompress(aorname(n) + strmid(list100(i),1))
   endfor
 openw, outlun200_64, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_200_64s.txt'),/GET_LUN
   for i = 0, count200_64 - 1 do begin
      printf, outlun200_64, list200_64(i)
      printf,58, strcompress(aorname(n) + strmid(list200_64(i),1))
   endfor

   openw, outlun100_32, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_100_32s.txt'),/GET_LUN
   for i = 0, count100_32 - 1 do begin
      printf, outlun100_32, list100_32(i)
      printf, 57, strcompress(aorname(n) + strmid(list100_32(i),1))
   endfor

   openw, outlun30, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_30s.txt'),/GET_LUN
   for i = 0, count30 - 1 do begin
      printf, outlun30, list30(i)
      printf, 56, strcompress(aorname(n) + strmid(list30(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 55, strcompress(aorname(n) + strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 54, strcompress(aorname(n) + strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 53, strcompress(aorname(n) + strmid(list2(i),1))
   endfor

   openw, outlun0pt6, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_0pt6s.txt'),/GET_LUN
   for i = 0, count0pt6 - 1 do begin
      printf, outlun0pt6, list0pt6(i)
      printf, 94, strcompress(aorname(n) + strmid(list0pt6(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 52, strcompress(aorname(n) + strmid(list0pt4(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/bcd/' + aorname(n)+ '/bcdlist_ch2_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 51, strcompress(aorname(n) + strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
 close, outlun200_64
   close, outlun100_32
   close, outlun30
   close, outlun6
   close, outlun2
   close, outlun0pt6
   close, outlun0pt4
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun12
   free_lun, outlun200_64
   free_lun, outlun100_32
   free_lun, outlun30
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt6
   free_lun, outlun0pt4
   free_lun, outlun0pt1


  
endfor

close, /all 

end
