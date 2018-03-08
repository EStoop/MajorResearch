pro pixel_dist, planetname
  
;;look at the first image in a set of AORs and make a histogram of the
;;pixel values in that BCD.

  colorarr = ['burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']  
  chname = '2'
  bsize = 1.0  ;;binsize for histograms

  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')

  for a = 0, n_elements(aorname) - 1 do begin
     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dirname+'bcdlist.txt')
     spawn, command
     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent

     ;;read in the first fits image
     fits_read, fitsname(4), im, h
     plothist, im, xhist, yhist,/NAN, /noprint, /noplot
     pa = plot(xhist, yhist, xtitle = 'mJy', ytitle = 'Number', color = colorarr(a), overplot = pa,/ylog, $
               title = planetname, xrange = [-100, 800])
  endfor


end
