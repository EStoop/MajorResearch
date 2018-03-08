function run_mopex, directory_name, mopex_script_env, namelist_ch1, namelist_ch2
;function run_mopex, directory_name, mopex_script_env, overlap_nl_ch1, overlap_nl_ch2,namelist_ch1, namelist_ch2

; example directory_name where the data live= '/Users/jkrick/IRAC/EGS/download'
; example mopex_script = '~/bin/mopex/mopex-script-env.csh'
; example namelist = 'mosaic_test_namelist.nl'

;make pertinent file lists for input to MOPEX
command1 = ' find '+directory_name+' -name "*.fits" > ' + directory_name + '/files_test.list'

command2 = 'grep ch1 < '+directory_name+'/files_test.list | grep -v pmask | grep _bcd.fits > '+directory_name+'/ch1_bcd_test.list'
command3 = 'grep ch1 < '+directory_name+'/files_test.list | grep _bunc.fits > '+directory_name+'/ch1_bunc_test.list'
command4 = 'grep ch1 < '+directory_name+'/files_test.list | grep _bdmsk.fits > '+directory_name+'/ch1_bimsk_test.list'
command5 = 'grep ch2 < '+directory_name+'/files_test.list | grep -v pmask | grep _bcd.fits > '+directory_name+'/ch2_bcd_test.list'
command6 = 'grep ch2 < '+directory_name+'/files_test.list | grep _bunc.fits > '+directory_name+'/ch2_bunc_test.list'
command7 = 'grep ch2 < '+directory_name+'/files_test.list | grep _bdmsk.fits > '+directory_name+'/ch2_bimsk_test.list'

;get rid of the first frame in all of the directories
command8 = 'cat ' + directory_name + '/ch1_bcd_test.list | grep -v 0000_0000 > ' + directory_name+'/ch1_bcd_good.list'
command9 = 'cat ' + directory_name + '/ch1_bunc_test.list | grep -v 0000_0000 > ' + directory_name+'/ch1_bunc_good.list'
command10 = 'cat ' + directory_name + '/ch1_bimsk_test.list | grep -v 0000_0000 > ' + directory_name+'/ch1_bimsk_good.list'
command11 = 'cat ' + directory_name + '/ch2_bcd_test.list | grep -v 0000_0000 > ' + directory_name+'/ch2_bcd_good.list'
command12 = 'cat ' + directory_name + '/ch2_bunc_test.list | grep -v 0000_0000 > ' + directory_name+'/ch2_bunc_good.list'
command13 = 'cat ' + directory_name + '/ch2_bimsk_test.list | grep -v 0000_0000 > ' + directory_name+'/ch2_bimsk_good.list'


;run these commands
commands = [command1, command2, command3, command4, command5, command6, command7, command8, command9, command10, command11, command12, command13]

for i = 0, n_elements(commands)-1 do spawn, commands(i)

;set environmental variables and run mopex
command14 = 'source ' + mopex_script_env

command15 = 'overlap.pl -n ' + namelist_ch1 + ' > overlap_logfile_ch1.txt'
command16 = 'overlap.pl -n ' + namelist_ch2 + ' > overlap_logfile_ch2.txt'
command17 = 'mosaic.pl -n ' + namelist_ch1 + ' > mosaic_logfile_ch1.txt'
command18 = 'mosaic.pl -n ' + namelist_ch2 + ' > mosaic_logfile_ch2.txt'

CD, directory_name

commands = 'pwd ; ' + command14 + ' ; ' + command17 + ' ; ' + command18
spawn, commands

return, 0
end
