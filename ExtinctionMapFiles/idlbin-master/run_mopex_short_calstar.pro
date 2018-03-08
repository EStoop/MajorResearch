function run_mopex_short_calstar, directory_name, mopex_script_env, namelist_ch1, namelist_ch2, framtime

full_name = directory_name  ;full directory path to calaor
;calaor= strmid(directory_name, 49)  ;cal aor name
directory_name = strmid(directory_name,0,37) ;path just to 016300 ;37

print,'full_name', full_name
print, 'directory_name', directory_name

; example directory_name where the data live= '/Users/jkrick/IRAC/EGS/download'
; example mopex_script = '~/bin/mopex/mopex-script-env.csh'
; example namelist = 'mosaic_test_namelist.nl'

;make pertinent file lists for input to MOPEX
command1 = ' find '+full_name+' -name "*.fits" > ' + directory_name + '/files_test.list'
print, 'command1', command1

command2 = 'grep IRAC.1. < '+directory_name+'/files_test.list | grep .bcd_fp.fits > '+directory_name+'/ch1_bcd_test.list'
command3 = 'grep IRAC.1. < '+directory_name+'/files_test.list | grep flux.fits > '+directory_name+'/ch1_bunc_test.list'
command4 = 'grep IRAC.1. < '+directory_name+'/files_test.list | grep dmask.fits > '+directory_name+'/ch1_bimsk_test.list'
command5 = 'grep IRAC.2. < '+directory_name+'/files_test.list | grep .bcd_fp.fits > '+directory_name+'/ch2_bcd_test.list'
command6 = 'grep IRAC.2. < '+directory_name+'/files_test.list | grep flux.fits > '+directory_name+'/ch2_bunc_test.list'
command7 = 'grep IRAC.2. < '+directory_name+'/files_test.list | grep dmask.fits > '+directory_name+'/ch2_bimsk_test.list'

;get rid of the first frame in all of the directories
command8 = 'cat ' + directory_name + '/ch1_bcd_test.list | grep -v 0000.0000  | grep -v flat > ' + directory_name+'/ch1_bcd_good.list'  ;| grep  002
command9 = 'cat ' + directory_name + '/ch1_bunc_test.list | grep -v 0000.0000    > ' + directory_name+'/ch1_bunc_good.list'  ;| grep  002
command10 = 'cat ' + directory_name + '/ch1_bimsk_test.list | grep -v 0000.0000  > ' + directory_name+'/ch1_bimsk_good.list' ;| grep  002
command11 = 'cat ' + directory_name + '/ch2_bcd_test.list | grep -v 0000.0000  | grep -v flat > ' + directory_name+'/ch2_bcd_good.list' ;| grep  004
command12 = 'cat ' + directory_name + '/ch2_bunc_test.list | grep -v 0000.0000 > ' + directory_name+'/ch2_bunc_good.list';| grep  004
command13 = 'cat ' + directory_name + '/ch2_bimsk_test.list | grep -v 0000.0000 > ' + directory_name+'/ch2_bimsk_good.list';| grep  004


;run these commands
commands = [command1, command2, command3, command4, command5, command6, command7, command8, command9, command10, command11, command12, command13]

for i = 0, n_elements(commands)-1 do spawn, commands(i)


;pick out just the ones with the exposure time we are interested in:
;readcol, directory_name + '/ch1_bcd_good.list', bcd_name, format= "A", /silent
;readcol, directory_name + '/ch1_bunc_good.list', bunc_name, format= "A", /silent
;readcol, directory_name + '/ch1_bimsk_good.list', bimsk_name, format= "A", /silent
;goodbcd = bcd_name
;goodbunc = bunc_name
;goodbimsk = bimsk_name
;m = 0
;


;read in the bcd fits files and multiply by the array location
;dependent flat field correction from SSC website
;fits_read, '/Users/jkrick/iwic/photcorr/ch1_photcorr_rj.fits', corrdata, corrheader


;for k = 0, n_elements(bcd_name) -1 do begin

;   bcdheader = headfits(bcd_name(k)) ;
;   frametime = fxpar(bcdheader, 'framtime')
;   chan = fxpar(bcdheader, 'chnlnum')
;   if frametime eq framtime and chan eq 1 then begin
;      fits_read,  bcd_name(k), data, header
;      data = data * corrdata
;      newfilename = strmid(bcd_name(k), 0, 52) + '/flat' + strmid(bcd_name(k), 53)
      ;print, 'newfilename', newfilename
      ;fits_write, newfilename, data, header

;      goodbcd(m) = newfilename
;      goodbunc(m) = bunc_name(k)
;      goodbimsk(m) = bimsk_name(k)
;      m = m + 1
;   endif
;endfor



;need to write out the files for mopex to read back in

;openw, outlunbcd, directory_name + '/ch1_bcd_good.list', /GET_LUN
;openw, outlunbunc, directory_name + '/ch1_bunc_good.list', /GET_LUN
;openw, outlunbimsk, directory_name + '/ch1_bimsk_good.list', /GET_LUN

;in this case I really only want to combine 7 of them
;for n = 0, m-1 do begin ;m-1
;   printf, outlunbcd, goodbcd(n)
;   printf, outlunbunc, goodbunc(n)
;   printf, outlunbimsk, goodbimsk(n)
;endfor

;close, outlunbcd
;free_lun, outlunbcd
;close, outlunbunc
;free_lun, outlunbunc
;close, outlunbimsk
;free_lun, outlunbimsk



;--------------------------------------
;ch2 

;pick out just the ones with the exposure time we are interested in:
;readcol, directory_name + '/ch2_bcd_good.list', bcd_name, format= "A", /silent
;readcol, directory_name + '/ch2_bunc_good.list', bunc_name, format= "A", /silent
;readcol, directory_name + '/ch2_bimsk_good.list', bimsk_name, format= "A", /silent
;goodbcd = bcd_name
;goodbunc = bunc_name
;goodbimsk = bimsk_name
;m = 0


;read in the bcd fits files and multiply by the array location
;dependent flat field correction from SSC website
;fits_read, '/Users/jkrick/iwic/photcorr/ch2_photcorr_rj.fits', corrdata, corrheader


;for k = 0, n_elements(bcd_name) -1 do begin
;   bcdheader = headfits(bcd_name(k)) ;
;   frametime = fxpar(bcdheader, 'framtime')
;   chan = fxpar(bcdheader, 'chnlnum')

;   if frametime eq framtime and chan eq 2 then begin
;      fits_read,  bcd_name(k), data, header
;      data = data * corrdata
;      newfilename = strmid(bcd_name(k), 0, 52) + '/flat' + strmid(bcd_name(k), 53)
;;      print, 'goodbcd(count)', newfilename
;;      fits_write, newfilename, data, header

;      goodbcd(m) = newfilename
;;      goodbunc(m) = bunc_name(k)
;      goodbimsk(m) = bimsk_name(k)
;      m = m + 1
;   endif
;endfor

;need to write out the files for mopex to read back in

;openw, outlunbcd, directory_name + '/ch2_bcd_good.list', /GET_LUN
;;openw, outlunbunc, directory_name + '/ch2_bunc_good.list', /GET_LUN
;openw, outlunbimsk, directory_name + '/ch2_bimsk_good.list', /GET_LUN
;;in this case I really only want to combine 7 of them
;for n = 0, m-1 do begin ;m - 1
;   printf, outlunbcd, goodbcd(n)
;;   printf, outlunbunc, goodbunc(n)
;   printf, outlunbimsk, goodbimsk(n)
;endfor

;close, outlunbcd
;free_lun, outlunbcd
;;close, outlunbunc
;;free_lun, outlunbunc
;close, outlunbimsk
;free_lun, outlunbimsk

;--------------------------------------

;set environmental variables and run mopex
command14 = 'source ' + mopex_script_env

command15 = 'overlap.pl -n ' + namelist_ch1 + ' > overlap_logfile_ch1.txt'
command16 = 'overlap.pl -n ' + namelist_ch2 + ' > overlap_logfile_ch2.txt'
command17 = 'mosaic.pl -n ' + namelist_ch1 + ' > mosaic_logfile_ch1.txt'
command18 = 'mosaic.pl -n ' + namelist_ch2 + ' > mosaic_logfile_ch2.txt'

;up_dir = strmid(directory_name, 0, 35)  ;don't want to be in AOR
;directory, cal files aren't there.
print, 'cd', directory_name
CD, directory_name;up_dir

commands = 'pwd ; ' + command14 + ' ; ' + command17 + ' ; ' + command18
print, 'starting mopex'
spawn, commands
print, 'done running mopex'
return, 0
end
