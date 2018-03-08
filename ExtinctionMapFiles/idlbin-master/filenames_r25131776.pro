pro filenames, dir_name

;dir_name = '/Users/jkrick/iwic/r25131776'

aorname = 'r25131776'
exptime = [0.1,0.4,2,6,12,100,200]

;need to create a rawlist; ch1, position 1, exptime
;readcol, dir_name + '/aor_list', aorname, format="A",/silent

;open file lists for all three positions together

openw,90, strcompress(dir_name+ '/ch1/raw/' + 'rawlist_ch1_200s.txt')
openw,89, strcompress(dir_name+ '/ch1/raw/' + 'rawlist_ch1_100s.txt')
openw,88, strcompress(dir_name+ '/ch1/raw/' + 'rawlist_ch1_12s.txt')
openw,87, strcompress(dir_name+ '/ch1/raw/' +  'rawlist_ch1_6s.txt')
openw,86, strcompress(dir_name+ '/ch1/raw/' +  'rawlist_ch1_2s.txt')
openw,85, strcompress(dir_name+ '/ch1/raw/' + 'rawlist_ch1_0pt4s.txt')
openw,84, strcompress(dir_name+ '/ch1/raw/' + 'rawlist_ch1_0pt1s.txt')
openw,83, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_200s.txt')
openw,82, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_100s.txt')
openw,81, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_12s.txt')
openw,80, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_6s.txt')
openw,79, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_2s.txt')
openw,78, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_0pt4s.txt')
openw,77, strcompress(dir_name+ '/ch1/bcd/' + 'bcdlist_ch1_0pt1s.txt')
openw,76, strcompress(dir_name+ '/ch2/raw/' + 'rawlist_ch2_200s.txt')
openw,75, strcompress(dir_name+ '/ch2/raw/' + 'rawlist_ch2_100s.txt')
openw,74, strcompress(dir_name+ '/ch2/raw/' + 'rawlist_ch2_12s.txt')
openw,73, strcompress(dir_name+ '/ch2/raw/' +  'rawlist_ch2_6s.txt')
openw,72, strcompress(dir_name+ '/ch2/raw/' +  'rawlist_ch2_2s.txt')
openw,71, strcompress(dir_name+ '/ch2/raw/' + 'rawlist_ch2_0pt4s.txt')
openw,70, strcompress(dir_name+ '/ch2/raw/' + 'rawlist_ch2_0pt1s.txt')
openw,69, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_200s.txt')
openw,68, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_100s.txt')
openw,67, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_12s.txt')
openw,66, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_6s.txt')
openw,65, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_2s.txt')
openw,64, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_0pt4s.txt')
openw,63, strcompress(dir_name+ '/ch2/bcd/' + 'bcdlist_ch2_0pt1s.txt')



;for each of the AORs (3)
for n = 0,  n_elements(aorname) - 1 do begin


;   print, 'working on aor', aorname(n)
   ;change to that AOR directory
;make lists by channel number for raw and bcd
   CD, strcompress(dir_name + '/ch1/raw/')
   command1 =  ' find . -name "*.fits" >  raw_list.txt'
   command2 = 'cat raw_list.txt | grep I1. > raw_ch1.txt'
   commands = [command1, command2]
   for i = 0, n_elements(commands) -1 do spawn, commands(i)

   CD, strcompress(dir_name + '/ch2/raw/')
   command1 =  ' find . -name "*.fits" >  raw_list.txt'
   command3 = 'cat raw_list.txt | grep I2. > raw_ch2.txt'
   commands = [command1, command3]
   for i = 0, n_elements(commands) -1 do spawn, commands(i)

;now bcd
 
   CD, strcompress(dir_name + '/ch1/bcd/')
   command1 =  ' find . -name "*bcd.fits" >  bcd_list.txt'
   command2 = 'cat bcd_list.txt | grep I1. > bcd_ch1.txt'
   commands = [command1, command2]
   for i = 0, n_elements(commands) -1 do spawn, commands(i)

   CD, strcompress(dir_name + '/ch2/bcd/')
   command1 =  ' find . -name "*bcd.fits" >  bcd_list.txt'
   command3 = 'cat bcd_list.txt | grep I2. > bcd_ch2.txt'
   commands = [command1, command3]
   for i = 0, n_elements(commands) -1 do spawn, commands(i)

;now need to make lists by exposure time
;which means reading in the headers of these fits files.  all of them
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;ch1
;start with the raws
   readcol, strcompress(dir_name + '/ch1/raw/raw_ch1.txt'), rawch1list, format="A", /silent
   count200 = 0
   list200 = strarr(2000)
   count100 = 0
   list100 = strarr(2000)
   count12 = 0
   list12 = strarr(2000)
   count6 = 0
   list6 = strarr(2000)
   count2 = 0
   list2 = strarr(2000)
   count0pt4 = 0
   list0pt4 = strarr(2000)
   count0pt1 = 0
   list0pt1 = strarr(2000)

   CD, strcompress(dir_name + '/ch1/raw/')

   for j = 0, n_elements(rawch1list) - 1 do begin
      rawheader = headfits(rawch1list(j)) ;
      fowlnum = fxpar(rawheader, 'A0614D00')
      waitper = fxpar(rawheader, 'A0615D00')

      if fowlnum eq 32 then begin
         list200(count200) = rawch1list(j)
         count200 = count200+1
      endif
      if fowlnum eq 16 and waitper eq 468 then begin
         list100(count100) = rawch1list(j)
         count100 = count100+1
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
      if fowlnum eq 1 and waitper eq 0 then begin
         list0pt4(count0pt4) = rawch1list(j)
         count0pt4 = count0pt4+1
      endif
      if fowlnum eq 2 and waitper eq 6 then begin
         list0pt1(count0pt1) = rawch1list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor


   openw, outlun200, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,90, strcompress(  strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 89, strcompress(  strmid(list100(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 88, strcompress(  strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 87, strcompress(  strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 86, strcompress(  strmid(list2(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 85, strcompress(  strmid(list0pt4(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/ch1/raw/'  + '/rawlist_ch1_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 84, strcompress(  strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
   close, outlun6
   close, outlun2
   close, outlun0pt4
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun12
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt4
   free_lun, outlun0pt1
   

;----------------------------------------------------------------------------------
;----------------------------------------------------------------------------------
;could fake the bcd's
;faster just to copy everything from above

;go on the bcd's
   readcol, strcompress(dir_name + '/ch1/bcd/bcd_ch1.txt'), bcdch1list, format="A", /silent
   count200 = 0
   list200 = strarr(2000)
   count100 = 0
   list100 = strarr(2000)
   count12 = 0
   list12 = strarr(2000)
   count6 = 0
   list6 = strarr(2000)
   count2 = 0
   list2 = strarr(2000)
   count0pt4 = 0
   list0pt4 = strarr(2000)
   count0pt1 = 0
   list0pt1 = strarr(2000)

   CD, strcompress(dir_name + '/ch1/bcd/' )

   for j = 0, n_elements(bcdch1list) - 1 do begin
      bcdheader = headfits(bcdch1list(j)) ;
      frametime = fxpar(bcdheader, 'framtime')

      if frametime eq 200 then begin
         list200(count200) = bcdch1list(j)
         count200 = count200+1
      endif
      if frametime eq 100 then begin
         list100(count100) = bcdch1list(j)
         count100 = count100+1
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
      if frametime eq 0.4 then begin
         list0pt4(count0pt4) = bcdch1list(j)
         count0pt4 = count0pt4+1
      endif
      if frametime eq 0.1 then begin
         list0pt1(count0pt1) = bcdch1list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor

  
 openw, outlun200, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,83, strcompress(  strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 82, strcompress(  strmid(list100(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 81, strcompress(  strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 80, strcompress(  strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 79, strcompress(  strmid(list2(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 78, strcompress(  strmid(list0pt4(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/ch1/bcd/'  + '/bcdlist_ch1_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 77, strcompress(  strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
   close, outlun6
   close, outlun2
   close, outlun0pt4
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun12
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt4
   free_lun, outlun0pt1


;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
;ch2
;start with the raws
   readcol, strcompress(dir_name + '/ch2/raw/raw_ch2.txt'), rawch2list, format="A", /silent
   count200 = 0
   list200 = strarr(2000)
   count100 = 0
   list100 = strarr(2000)
   count12 = 0
   list12 = strarr(2000)
   count6 = 0
   list6 = strarr(2000)
   count2 = 0
   list2 = strarr(2000)
   count0pt4 = 0
   list0pt4 = strarr(2000)
   count0pt1 = 0
   list0pt1 = strarr(2000)

   CD, strcompress(dir_name + '/ch2/raw/' )

   for j = 0, n_elements(rawch2list) - 1 do begin
      rawheader = headfits(rawch2list(j)) ;
      fowlnum = fxpar(rawheader, 'A0614D00')
      waitper = fxpar(rawheader, 'A0615D00')

      if fowlnum eq 32 then begin
         list200(count200) = rawch2list(j)
         count200 = count200+1
      endif
      if fowlnum eq 16 and waitper eq 468 then begin
         list100(count100) = rawch2list(j)
         count100 = count100+1
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
      if fowlnum eq 1 and waitper eq 0 then begin
         list0pt4(count0pt4) = rawch2list(j)
         count0pt4 = count0pt4+1
      endif
      if fowlnum eq 2 and waitper eq 6 then begin
         list0pt1(count0pt1) = rawch2list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor


 openw, outlun200, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,76, strcompress(  strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 75, strcompress(  strmid(list100(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 74, strcompress(  strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 73, strcompress(  strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 72, strcompress(  strmid(list2(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 71, strcompress(  strmid(list0pt4(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/ch2/raw/'  + '/rawlist_ch2_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 70, strcompress(  strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
   close, outlun6
   close, outlun2
   close, outlun0pt4
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun12
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt4
   free_lun, outlun0pt1
   

;----------------------------------------------------------------------------------
;----------------------------------------------------------------------------------
;could fake the bcd's
;faster just to copy everything from above

;go on the bcd's
   readcol, strcompress(dir_name + '/ch2/bcd/bcd_ch2.txt'), bcdch2list, format="A", /silent
   count200 = 0
   list200 = strarr(2000)
   count100 = 0
   list100 = strarr(2000)
   count12 = 0
   list12 = strarr(2000)
   count6 = 0
   list6 = strarr(2000)
   count2 = 0
   list2 = strarr(2000)
   count0pt4 = 0
   list0pt4 = strarr(2000)
   count0pt1 = 0
   list0pt1 = strarr(2000)

   CD, strcompress(dir_name + '/ch2/bcd/' )

   for j = 0, n_elements(bcdch2list) - 1 do begin
      bcdheader = headfits(bcdch2list(j)) ;
      frametime = fxpar(bcdheader, 'framtime')

      if frametime eq 200 then begin
         list200(count200) = bcdch2list(j)
         count200 = count200+1
      endif
      if frametime eq 100 then begin
         list100(count100) = bcdch2list(j)
         count100 = count100+1
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
      if frametime eq 0.4 then begin
         list0pt4(count0pt4) = bcdch2list(j)
         count0pt4 = count0pt4+1
      endif
      if frametime eq 0.1 then begin
         list0pt1(count0pt1) = bcdch2list(j)
         count0pt1 = count0pt1+1
      endif

      
   endfor


openw, outlun200, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_200s.txt'),/GET_LUN
   for i = 0, count200 - 1 do begin
      printf, outlun200, list200(i)
      printf,69, strcompress(  strmid(list200(i),1))
   endfor

   openw, outlun100, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_100s.txt'),/GET_LUN
   for i = 0, count100 - 1 do begin
      printf, outlun100, list100(i)
      printf, 68, strcompress(  strmid(list100(i),1))
   endfor

   openw, outlun12, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_12s.txt'),/GET_LUN
   for i = 0, count12 - 1 do begin
      printf, outlun12, list12(i)
      printf, 67, strcompress(  strmid(list12(i),1))
   endfor

   openw, outlun6, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_6s.txt'),/GET_LUN
   for i = 0, count6 - 1 do begin
      printf, outlun6, list6(i)
      printf, 66, strcompress(  strmid(list6(i),1))
   endfor

   openw, outlun2, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_2s.txt'),/GET_LUN
   for i = 0, count2 - 1 do begin
      printf, outlun2, list2(i)
      printf, 65, strcompress(  strmid(list2(i),1))
   endfor

   openw, outlun0pt4, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_0pt4s.txt'),/GET_LUN
   for i = 0, count0pt4 - 1 do begin
      printf, outlun0pt4, list0pt4(i)
      printf, 64, strcompress(  strmid(list0pt4(i),1))
   endfor
      
   openw, outlun0pt1, strcompress(dir_name+ '/ch2/bcd/'  + '/bcdlist_ch2_0pt1s.txt'),/GET_LUN
   for i = 0, count0pt1 - 1 do begin
      printf, outlun0pt1, list0pt1(i)
      printf, 63, strcompress(  strmid(list0pt1(i),1))
   endfor

   close, outlun200
   close, outlun100
   close, outlun12
   close, outlun6
   close, outlun2
   close, outlun0pt4
   close, outlun0pt1
   
   free_lun, outlun200
   free_lun, outlun100
   free_lun, outlun12
   free_lun, outlun6
   free_lun, outlun2
   free_lun, outlun0pt4
   free_lun, outlun0pt1




  
endfor

close, /all 

end
