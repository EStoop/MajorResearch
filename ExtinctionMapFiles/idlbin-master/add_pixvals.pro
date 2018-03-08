pro add_pixvals, planetname, chname, apradius
  
;add pixel values in a 7x7x63 box around the central pixel to planethash


  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all) ;
  print, 'restoring ', savefilename
  restore, savefilename
  
  startaor = 0
  stopaor = n_elements(aorname) - 1
  for a = startaor,stopaor do begin

     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*bcd.fits' > "+dirname+'bcdlist.txt')
     spawn, command
     
     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     print, 'nfits', n_elements(fitsname), '  ', aorname(a)
     piarr = findgen(7,7,63*n_elements(fitsname)) ; ignore the first frame
;     help, piarr
     for i =0L,  n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
        
                                ;read in the files
        fits_read, fitsname(i), im, h
;        print, 'i', i,i*63
        piarr[i*63]= im[12:18, 12:18,1:*] 
     endfor                     ; for each fits file
     
     ;check for nan's
     bad = where(finite(piarr) lt 1, nbad)
     if nbad gt 0 then print, 'got a nan', nbad
     ;alter the planethash
     planethash[aorname(a),'pixvals'] = piarr

  endfor                        ; for each AOR
  save, planethash, filename=savefilename

  
end
