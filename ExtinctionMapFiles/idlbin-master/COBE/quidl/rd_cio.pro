pro rd_cio,file,product,field,data=data_list,pixel=pixel,startpix=startpix,$
     indfile=indfile

if not keyword_set(indfile) then begin
  case !VERSION.OS of
   'windows': strsep='\'
   'vms': strsep=']'
    else: strsep='/'
  endcase
  dirspec = strmid(file,0,rstrpos(file,strsep)+1)
  face = strpos(strupcase(file),'.FITS') 
  face = strlowcase(strmid(file,face-5,5))
  iname ='dirbe_cio_index_p3b_'+face+'.fits'
  indexname = dirspec + iname
  res = is_fits(indexname)
  if not res then indexname = dirspec + strsep + strupcase(iname)
endif else indexname = indfile
      
if not is_fits(indexname) then begin
     print,'Can not find index file: ',indexname
     read,indexname,prompt='Please enter index file name: '
endif
all_data = readfits(indexname,bin_header,/ext)
pixel = tbget(bin_header,all_data,1)
lastrow = tbget(bin_header,all_data,2)
all_data = 0
dum = temporary(all_data)

if n_elements(startpix) eq 0 then read,startpix,prompt='Enter pixel num:'
start_point = where(pixel eq startpix)
if start_point(0) eq -1 then begin
    message,/info,'You entered a pixel that was not observed this day'
    return
endif

if start_point(0) ne 0 then startrow = lastrow(start_point(0) -1) $
    else startrow = 0
stoprow = lastrow(start_point(0)) - 1
numrow = stoprow - startrow + 1

all_data = readfits(file,bin_header,/ext,startrow=startrow,numrow=numrow)

xfields = fxpar(bin_header,'TTYPE*')
unitlist = fxpar(bin_header,'TUNIT*')
if field eq '?' then begin
   print,'This file contains the following fields:'
   for i = 0, n_elements(xfields)-1 do print,xfields(i)
   field=''
   read,field,prompt='Enter Field Name as shown: '
endif
while strlen(field) lt 8 do field = field + ' '
choice = where(STRUPCASE(field) eq STRUPCASE(xfields))
unitlist = unitlist(choice(0))

data_list = tbget(bin_header,all_data,choice(0)+1)
if strupcase(field) ne 'PIXEL_NO' then $
     pixel = tbget(bin_header,all_data,1) else pixel = data_list
all_data = 0
dum = temporary(all_data)

end

