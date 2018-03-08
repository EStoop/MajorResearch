pro diff_darks_test
  
;
base_dir = '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/calibration/'

       filenm_a = strcompress(base_dir + 'r48687616/ch2/cal/SPITZER_I2_48998144_0000_1_C9737980_sdark.fits',/remove_all)
      ; header_a = headfits(filenm_a) 
       fits_read, filenm_a, data_a, header_a
       e_a = fxpar(header_a, 'FRAMTIME')
       ch_a = fxpar(header_a, 'CHNLNUM')

       filenm_b = strcompress(base_dir + 'r48688384/ch2/cal/SPITZER_I2_48987392_0000_1_C9735186_sdark.fits',/remove_all) 
       ;header_b = headfits(filenm_b) 
       fits_read, filenm_b, data_b, header_b
       e_b = fxpar(header_b, 'FRAMTIME')
       ch_b = fxpar(header_b, 'CHNLNUM')

       ;print, 'differenconfg', e_a, e_b, ch_a, ch_b
       
       diff = data_a / data_b
       ;print,   strcompress(base_dir+campaign_name +'_diff_' + string(ch_a) +'_'+ string(e_a) + '.fits',/remove_all)
       fits_write, strcompress(base_dir+'div.fits'), diff, header_a
       
 
end
