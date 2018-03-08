pro rd_dcaf,file,product,field,data=data_list,pixel=pixel,startpix=startpix,$
     stoppix=stoppix,indfile=indfile

if not keyword_set(indfile) then begin
  case !VERSION.OS of
   'windows': strsep='\'
   'vms': strsep=']'
    else: strsep='/'
  endcase
  dirspec = strmid(file,0,rstrpos(file,strsep)+1)
  face = strpos(strupcase(file),'FACE')
  face = strlowcase(strmid(file,face,5))
  case 1 of
   strpos(product,'DCAF') ne -1: iname ='dirbe_dcafindex_p3b_'+face+'.fits'
   strpos(product,'DSZA') ne -1: iname ='dirbe_dszaindex_'+face+'.fits'
  endcase
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
range1 = pixel(0)
range2 = pixel(n_elements(pixel)-1)
lastrow = tbget(bin_header,all_data,2)
all_data = 0
dum = temporary(all_data)

if n_elements(startpix) eq 0 then read,startpix,prompt='Enter start pixel'+$
     ' pixel range ='+strtrim(range1)+'-'+strtrim(range2)+': '
if startpix lt range1 or startpix gt range2 then begin
   print,'You entered an out of range pixel'
   read,startpix,prompt='Enter start pixel: pixel range ='+strtrim(range1)+$
     strtrim(range2)
endif
if n_elements(stoppix) eq 0 then read,stoppix,prompt='Enter stop pixel'+$
     ' pixel range ='+strtrim(range1)+'-'+strtrim(range2)+': '
if stoppix lt range1 or stoppix gt range2 then begin
   print,'You entered an out of range pixel'
   read,stoppix,prompt='Enter stop pixel: pixel range ='+strtrim(range1)+$
     strtrim(range2)
endif

start_point = where(pixel eq startpix)
stop_point = where(pixel eq stoppix)
if start_point(0) ne 0 then startrow = lastrow(start_point(0) -1) $
    else startrow = 0
stoprow = lastrow(stop_point(0)) - 1
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

if strupcase(strtrim(field,2)) eq 'STDDEV' then begin
         sz=size(data_list)
         dbsig=fltarr(sz(1),sz(2))
         print,'unpacking standard deviations'
         For ib=0,9 do begin
                     itemp = data_list(ib,*)
                     dim=SIZE(itemp)
                     CASE ib OF
                         4: n=3
                         5: n=3
                         6: n=2
                         7: n=2
                         8: n=0
                         9: n=1
                         else: n=4
                     ENDCASE
                     for k=0l,dim(2)-1 do $ 
                       CASE itemp(k) OF
                         0: dbsig(ib,k) = 10.0^(-n) 
                         255: dbsig(ib,k) = 10.0^(-n+4) 
                         else: dbsig(ib,k) = $
                                  10.0^((itemp(k) - 0.5)/63.5 - n)
                       ENDCASE
         endfor
         data_list=dbsig
         dbsig=0
endif

if strupcase(strtrim(field,2)) eq 'STOKESSD' then begin
         sz=size(data_list)
         dbsig=fltarr(sz(1),sz(2))
         n=4
         print,'unpacking Stokes standard deviations'
         For ib=0,5 do begin
                     itemp = data_list(ib,*)
                     dim=SIZE(itemp)
                     for k=0l,dim(2)-1 do $ 
                       CASE itemp(k) OF
                         0: dbsig(ib,k) = 10.0^(-n) 
                         255: dbsig(ib,k) = 10.0^(-n+4) 
                         else: dbsig(ib,k) = $
                                  10.0^((itemp(k) - 0.5)/63.5 - n)
                       ENDCASE
         endfor
         data_list=dbsig
         dbsig=0
endif


end

