pro make_raw_cube, rfiles, ofile

	n = n_elements(rfiles)
	rstack = fltarr(32, 32, 64, n)
	mstack = fltarr(n*64)
	for i = 0, n-1 do begin
		print_pct, i, n
		im = readfits(rfiles[i], h, /SILENT)
		bs = sxpar(h, 'A0741D00')
		fn = sxpar(h, 'A0614D00')
		wt = sxpar(h, 'A0615D00')
		readmod = sxpar(h, 'A0617D00')
		ch = sxpar(h, 'A0610D00')
; Flip sign		
		nim = 65535 - uint(im)

		nim = float(nim)
		nim = nim - 0.5 * (1. - 2.^(-bs))
		ptr = where(nim gt 45000., count)
		if (count gt 0) then nim[ptr] = nim[ptr] - 65535.
;		cube = cube * 2.^bs / fn
; turn into cube
		cube = float(reform(nim, 32, 32, 64))
		rstack[*, *, *, i] = cube
		for k = 0, 63 do mstack[i*63+k] = max(cube[*, *, k], /NAN)
	endfor
	save, FILENAME=ofile, rstack, mstack
return
end

pro deadband_raw_stack_driver
	dir = '/iracdata/flight/IWIC/IRAC026000/raw/'
	rfiles1 = dir + '0040837888/*.mipl.fits'
	rfiles2 = dir + '0040838144/*.mipl.fits'
	
	files1 = file_search(rfiles1)
	make_raw_cube, files1, 'r40837888_raw_cube.sav'
	files1 = file_search(rfiles2)
	make_raw_cube, files1, 'r40838144_raw_cube.sav'	
return
end

pro make_raw_maxdn_plots
return
end

pro make_simple_savefiles
	restore, 'prepatch_hd158460_ch2.sav'
	tmin = min(t2)
	ta = t2
	t0a = (ta - tmin) / 3600.
	xa = x32
	ya = y32
	sxa = ex32
	sya = ey32
	fa = f2[*,9] * 1.1226620D
	efa = ef2[*, 9] * 1.1226620D
	sfa = sf2
	
	restore, 'postpatch_hd158460_ch2.sav'
	tb = t2
	t0b = (tb - tmin) / 3600.
	xb = x32
	yb = y32
	sxb = ex32
	syb = ey32
	fb = f2[*,9] * 1.1226620D
	efb = ef2[*, 9] * 1.1226620D
	sfb = sf2
	
	ofile = 'deadband_photo.sav'
	save, FILENAME=ofile, ta, t0a, xa, ya, sxa, sya, fa, efa, sfa, $
	        tb, t0b, xb, yb, sxb, syb, fb, efb, sfb

return
end

pro compare_runs
	restore, 'deadband_photo.sav'
	tek_color
	!P.BACKGROUND=1
	!P.COLOR=0
	!P.PSYM=-3
	!P.THICK=1.25
	!P.CHARSIZE=1.5
	!P.SYMSIZE=1.75
	
; Use 100 for bottom and top borders, 100 for left and right, 300 for each plot 	
	window, /free, retain=2, xsize=1000, ysize=1100
	xt = 'Time from AOR start (hours)'
	yt1 = '(flux - median[flux]) / median[flux]'
	yt2 = 'X centroid - median(centroid)'
	yt3 = 'Y  centroid - median(centroid)'	
	tt = 'Deadband Test'
	tb1 = t0b - min(t0b)
;	plot, t0a, (fa-median(fa)) / median(fa), xtitle=xt, ytitle=yt1, xtick_get=xticks, $
;	      position=[100, 100, 900, 400], /DEVICE, PSYM=3, xrange=[0, 3.], yrange=[-0.015, 0.015],$
;	      /NOPLOT
	xtn = [' ', '-0.010', '-0.005', '0.000', '0.005', '0.010', ' ']
	plot, t0a, (fa-median(fa)) / median(fa), xtitle=xt, ytitle=yt1, xtick_get=xticks, $
	      position=[100, 100, 900, 400], /DEVICE, PSYM=3, xrange=[0, 3.], yrange=[-0.005, 0.005],$
	      /NODATA
	xyouts, 0.1, 0.003, 'Smoothed by 51 samples', CHARSIZE=2
	legend, ['Pre-reduction', 'Post-reduction'], COLORS=[0, 3], PSYM=[1, 1], /RIGHT, /TOP
;	oplot, tb1, (fb-median(fb)) / median(fb), color=3, psym=3
	favec = smooth((fa-median(fa)) / median(fa), 51, /EDGE_TRUNCATE)
	index = lindgen(488) * 51 + 25	
	oplot, t0a[index], favec[index], PSYM=1
	fbvec = smooth((fb-median(fb)) / median(fb), 51, /EDGE_TRUNCATE)
	oplot, t0a[index], fbvec[index], PSYM=1, COLOR=3
	blank = replicate(' ', n_elements(xticks))
;	plot, t0a, xa-median(xa), ytitle=yt2, xtitle='', xtickname=blank, $
;	      position=[100, 400, 900, 700], /DEVICE, /NOERASE, PSYM=3, $
;	      xrange=[0, 3.], yrange=[-0.1, 0.1], /noplot
	xtn = [' ', '-0.05', '0.00', '0.05', ' ']
	plot, t0a, xa-median(xa), ytitle=yt2, xtitle='', xtickname=blank, $
	      position=[100, 400, 900, 700], /DEVICE, /NOERASE, PSYM=3, $
	      xrange=[0, 3.], yrange=[-0.1, 0.1], ytickname=xtn 
	oplot, tb1, xb-median(xb), color=3, psym=3
;	plot, t0a, ya-median(ya), ytitle=yt3, xtitle='', xtickname=blank, $
;	      position=[100, 700, 900, 1000], /DEVICE, /NOERASE, PSYM=3, $
;	      xrange=[0, 3.], yrange=[-0.1, 0.1], title=tt, /noplot
	plot, t0a, ya-median(ya), ytitle=yt3, xtitle='', xtickname=blank, $
	      position=[100, 700, 900, 1000], /DEVICE, /NOERASE, PSYM=3, $
	      xrange=[0, 3.], yrange=[-0.1, 0.1], title=tt, ytickname=xtn 
	oplot, tb1, yb-median(yb), color=3, psym=3

	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'deadband_photo_raw.png', goo, r, g, b
	
    set_plot, 'PS'
    tek_color
    !P.BACKGROUND=1
    !P.COLOR=0
    !P.PSYM=-3

	!P.THICK=6.
    !P.CHARSIZE=2
    !P.CHARTHICK=1.25
    !P.SYMSIZE=6.
    !X.THICK=1.5
    !Y.THICK=1.5
    et = 5.0
    set_plot,'ps'
    !P.FONT = 0
    cg = 0
    mu = string(181B)
    sig = '!9' + string(115B) + '!6'
    xtime = '!9' + string(180B) + '!6'
    !P.FONT = 0
    fname='relphoto.eps'
    device,file=fname,xs=9,ys=9,/inches,/portrait,xo=.5,yo=.5,$
               /helvetica, /isolatin1,/encapsulated, /color

	plot, t0a, (fa-median(fa)) / median(fa), xtitle=xt, ytitle=yt1, xtick_get=xticks, $
	      PSYM=3, xrange=[0, 3.], yrange=[-0.005, 0.005],$
	      /NODATA
	xyouts, 0.1, 0.003, 'Smoothed by 101 samples', CHARSIZE=1.75
	legend, ['Pre-reduction', 'Post-reduction'], COLORS=[0, 3], PSYM=[-1, -1], /RIGHT, /TOP, $
	CHARSIZE=1.5
;	oplot, tb1, (fb-median(fb)) / median(fb), color=3, psym=3
	favec = smooth((fa-median(fa)) / median(fa), 101, /EDGE_TRUNCATE)
	index = lindgen(246) * 101 + 50	
	oplot, t0a[index], favec[index], PSYM=-1
	fbvec = smooth((fb-median(fb)) / median(fb), 101, /EDGE_TRUNCATE)
	oplot, t0a[index], fbvec[index], PSYM=-1, COLOR=3
               
        device, /CLOSE
		set_plot, 'X'
		!P.FONT=-1   
		stop
return
end

pro write_deadband_photo
	fstr = '(2(F14.3, 2(2x, F9.5, 2x, F7.5), 2(2x, F9.7), 2x))'
;##########.###  ##.#####  #.#####  ##.#####  #.#####  #.########  #.#######
	restore, 'deadband_photo.sav'
	openw, lun, 'deadband_photo.txt', /GET_LUN
	n = n_elements(ta)
	for i = 0, n-1 do $
	         printf, FORMAT=fstr, lun, ta[i], xa[i], sxa[i], ya[i], sya[i], fa[i], efa[i], $
	                                       tb[i], xb[i], sxb[i], yb[i], syb[i], fb[i], efb[i]
	free_lun, lun
return
end

pro fit_pp_driver
openw, lun, 'deadband.log', /GET_LUN

	tek_color
	!P.BACKGROUND=1
	!P.COLOR=0
	!P.PSYM=-3
	!P.THICK=1.25
	!P.CHARSIZE=1.5
	!P.SYMSIZE=1.75
	
; From previous fit, 
	asig = 0.0026556005
	bsig = 0.0025597424

	restore, 'deadband_photo.sav'
	np = 6
	fname = 'fpa1_xfunc2'
	
	window, 15, retain=2, xsize=500, ysize=500
	window, 16, retain=2, xsize=900, ysize=600	
	
; Remove pointing outliers
	axhi = 22.95 & ayhi = 231.37
	agptr = where(xa lt axhi and ya lt ayhi)
	aflag = xa * 0
	aflag[agptr] = 1
	goo = where(aflag eq 0, goocount)
	printf, lun, 'Number of centroid outliers removed in pre-data is ' + strn(goocount)
	
	bxlo = 22.981593 & bxhi = 23.120959 & bylo = 231.09644 & byhi = 231.24318
	bgptr = where(xb gt bxlo and xb lt bxhi and yb gt bylo and yb lt byhi)
	bflag = xb * 0
	bflag[bgptr] = 1
	goo = where(bflag eq 0, goocount)
	printf, lun, 'Number of centroid outliers removed in post-data is ' + strn(goocount)	
	bflo = 1.1773940

	xt = 'X (pixels)'
	yt = 'Y (pixels)'
	xt1 = 'Time (hours)'
	yt1 = 'Flux (Jy)'
	tta = 'r40837888 ch2 pre-patch'
	ttb = 'r40838144 ch2 post-patch'

	xra = [22.7, 23.9]
	yra = [231.0, 231.7]
	wset, 15
	plot, xa, ya, PSYM=1, xtitle=xt, ytitle=yt, title=tta, xrange=xra, $
	      yrange=yra
	plots, [axhi, axhi], yra, PSYM=-3, color=2
	plots, xra, [ayhi, ayhi], PSYM=-3, color=2
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'pre_xy_scatter.png', goo, r, g, b
	
	wset, 16
	plot, t0a, fa, PSYM=4, xtitle=xt1, ytitle=yt1, title=tta, yrange=[1.14, 1.33], $
	      xrange=[0., 3.]
	sfptra = where(sfa eq 0 and aflag eq 1, goocount)
	printf, lun, 'Dropped ' + strn(goocount) + ' samples as first subframe in pre-data'
	oplot, t0a[sfptra], fa[sfptra], PSYM=4, COLOR=4
	bptr = where(aflag eq 0)
	oplot, t0a[bptr], fa[bptr], PSYM=4, COLOR=2
	legend, ['Good', 'First subframe', 'Bad centroid'], PSYM=4, COLORS=[0, 4, 2], /RIGHT
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'pre_phot_time.png', goo, r, g, b	

	aflag[sfptra] = 0
	aptr = where(aflag eq 1)
	tta1 = tta + ' Trimmed'
	plot, t0a[aptr], fa[aptr], PSYM=4, xtitle=xt1, ytitle=yt1, title=tta1, yrange=[1.16, 1.23], $
	      xrange=[0., 3.]
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'pre_phot_time_trimmed.png', goo, r, g, b	
	
; Drop first three minutes of data and most significant outliers by eye
	ptr = where(t0a lt 0.05 and aflag eq 1, goocount)
	printf, lun, 'Dropped ' + strn(goocount) + ' samples at beginning of light curve in pre-data'
	aflag[ptr] = 0
	agptr = where(aflag eq 1)
	ma = median(fa[agptr])
	apct = (fa - ma) / ma
	abptr = where(abs(apct) gt 0.009 and aflag eq 1, goocount)
	printf, lun, 'Dropped ' + strn(goocount) + ' samples > 0.9% from median in pre-data'
	aflag[abptr] = 0
	
; Initial guess	
	agptr = where(aflag eq 1)
	pa0 = [median(fa[agptr]), fltarr(np-1)]
	amx = median(xa[agptr])
	amy = median(ya[agptr])	
	adx = xa[agptr] - amx
	ady = ya[agptr] - amy
	
;	afargs = {FLUX:fa[agptr], DX:adx, DY:ady, T:t0a[agptr], ERR:efa[agptr]}
	afargs = {FLUX:fa[agptr], DX:adx, DY:ady, T:t0a[agptr], ERR:efa[agptr]*0.+asig}
	pa = mpfit(fname, pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof)
	achi = achi / adof
	
	fpa1_model2, pa, xa-amx, ya-amy, t0a, amodel
	ct0a = t0a[agptr]
	cfa = fa[agptr] / amodel[agptr]
	capct = (cfa - ma) / ma
	
	xrb = [22.6, 23.2]
	yrb = [230.5, 231.4]
	wset, 15
	plot, xb, yb, PSYM=1, xtitle=xt, ytitle=yt, title=ttb, xrange=xrb, $
	      yrange=yrb
	plots, [bxhi, bxhi], yrb, PSYM=-3, color=2
	plots, xrb, [byhi, byhi], PSYM=-3, color=2
	plots, [bxlo, bxlo], yrb, PSYM=-3, color=2
	plots, xrb, [bylo, bylo], PSYM=-3, color=2	
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'post_xy_scatter.png', goo, r, g, b
	
	wset, 16
	plot, t0b, fb, PSYM=4, xtitle=xt1, ytitle=yt1, title=ttb, yrange=[1.15, 1.35], $
	      xrange=[2.9, 6.]
	sfptrb = where(sfb eq 0 and bflag eq 1, goocount)
	printf, lun, 'Dropped ' + strn(goocount) + ' samples as first subframe in post-data'	
	oplot, t0b[sfptrb], fb[sfptrb], PSYM=4, COLOR=4
	bptr = where(bflag eq 0)
	oplot, t0b[bptr], fb[bptr], PSYM=4, COLOR=2
	legend, ['Good', 'First subframe', 'Bad centroid'], PSYM=4, COLORS=[0, 4, 2], /RIGHT
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'post_phot_time.png', goo, r, g, b	

	bflag[sfptrb] = 0
	bptr = where(bflag eq 1)
	ttb1 = ttb + ' Trimmed'
	plot, t0b[bptr], fb[bptr], PSYM=4, xtitle=xt1, ytitle=yt1, title=ttb1, yrange=[1.17, 1.24], $
	      xrange=[2.9, 6.]
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'post_phot_time_trimmed.png', goo, r, g, b	
	
; Drop first three minutes of data and most significant outliers by eye
	ptr = where(t0b lt (min(t0b) + 3./60.) and bflag eq 1, goocount)	
	printf, lun, 'Dropped ' + strn(goocount) + ' samples at beginning of light curve in post-data'	
	bflag[ptr] = 0
	bgptr = where(bflag eq 1)
	mb = median(fb[bgptr])
	bpct = (fb - mb) / mb
	bbptr = where(abs(bpct) gt 0.007 and bflag eq 1, goocount)
	printf, lun, 'Dropped ' + strn(goocount) + ' samples > 0.9% from median in pre-data'	
	bflag[bbptr] = 0	
	
; Initial guess	
	bgptr = where(bflag eq 1)
	pb0 = [median(fb[bgptr]), fltarr(np-1)]
	bmx = median(xb[bgptr])
	bmy = median(yb[bgptr])	
	bdx = xb[bgptr] - bmx
	bdy = yb[bgptr] - bmy
	
	bfargs = {FLUX:fb[bgptr], DX:bdx, DY:bdy, T:t0b[bgptr], ERR:efb[bgptr]}
	bfargs = {FLUX:fb[bgptr], DX:bdx, DY:bdy, T:t0b[bgptr], ERR:efb[bgptr]*0.+bsig}	
	pb = mpfit(fname, pb0, FUNCTARGS=bfargs, PERROR=spb, BESTNORM=bchi, DOF=bdof)
	bchi = bchi / bdof
	
	fpa1_model2, pb, xb-bmx, yb-bmy, t0b, bmodel
	ct0b = t0b[bgptr]
	cfb = fb[bgptr] / bmodel[bgptr]
	cbpct = (cfb - mb) / mb	
	
	acfmom = moment(cfa)
	bcfmom = moment(cfb)
	amom = moment(capct)
	bmom = moment(cbpct)
		
	dmin = -0.007
	dmax = 0.007
	bs = 0.0001
	ad = histogram(capct, MIN=dmin, MAX=dmax, BINSIZE=bs, LOCATION=al)
	bd = histogram(cbpct, MIN=dmin, MAX=dmax, BINSIZE=bs, LOCATION=bl)
	plot, al + bs/2., ad, psym=10, xrange=[-0.0075, 0.0075], xtitle='(flux-median)/median', $
	      ytitle='N', title='Flux variation after correction', yrange=[0, 500],THICK=2.5, $
	      CHARSIZE=2.
	oplot, bl + bs/2., bd, psym=10, color=3, THICK=2.5
	legend, ['Pre', 'Post'], PSYM=[-3, -3], COLORS=[0, 3]
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'flux_var_correct.png', goo, r, g, b
	
	set_ps_values, 'flux_var_correct.eps', et, cg, mu, sig, xtime, /landscape
	plot, al + bs/2., ad, psym=10, xrange=[-0.0075, 0.0075], xtitle='(flux-median)/median', $
	      ytitle='N', title='Flux variation after correction', yrange=[0, 500],THICK=3.
	oplot, bl + bs/2., bd, psym=10, color=3, THICK=3.,SYMSIZE=3.
	legend, ['Pre', 'Post'], PSYM=[-3, -3], COLORS=[0, 3]	
    end_ps

	
	avg_stack, cfa[0:23999], adb, ars, ads, amw, asw
	avg_stack, cfb[0:23999], bdb, brs, bds, bmw, bsw
	
	plot_noise_vs_bin, adb, asw, bdb, bsw
	
	wset, 15
	tek_color
	!P.BACKGROUND=1
	!P.COLOR=0
	!P.PSYM=-3
	!P.THICK=1.25
	!P.CHARSIZE=1.5
	!P.SYMSIZE=1.75
	plot, t0a, amodel, xtitle='!3Time (hours)', ytitle='!3C(!7D!Xx, !7D!Xy)', title='!3Pre',$
	      yrange=[0.995, 1.005]
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'amodel.png', goo, r, g, b	
	
	plot, t0b-min(t0b), bmodel, xtitle='!3Time (hours)', ytitle='!3C(!7D!Xx, !7D!Xy)', $
	      title='!3Post', yrange=[0.998, 1.002]
	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'bmodel.png', goo, r, g, b	

; #.###### #.######  #.###### #.###### 
	printf, lun
	
; Try periodigrams to find amplitude and period of pixel-phase effect and pointing wobble
	bxnu = 1 & bxlnp = bxnu & bynu = bxnu & bylnp = bxnu
	bfnu = 1 & bflnp = bxnu & cbfnu = 1 & cbflnp = bxnu
	bxlnp_coef = lnp_test(t0b[bgptr], bdx, /double, wk1=bxnu, wk2=bxlnp)
	bylnp_coef = lnp_test(t0b[bgptr], bdy, /double, wk1=bynu, wk2=bylnp)
	bflnp_coef = lnp_test(t0b[bgptr],(fb[bgptr]-mb)/mb, /double, wk1=bfnu, wk2=bflnp)
	cbflnp_coef = lnp_test(t0b[bgptr],(cfb-mb)/mb, /double, wk1=cbfnu, wk2=cbflnp)

	
	bxlnp_ptr = where(bxlnp eq max(bxlnp))
	bxnu_max = bxnu[bxlnp_ptr]
	bylnp_ptr = where(bylnp eq max(bylnp))
	bynu_max = bynu[bylnp_ptr]	

	print, 'New period = ' + strn(2. / (bxnu_max + bynu_max))

	axnu = 1 & axlnp = axnu & aynu = axnu & aylnp = axnu
	afnu = 1 & aflnp = axnu & cafnu = 1 & caflnp = axnu	
	axlnp_coef = lnp_test(t0a[agptr], adx, /double, wk1=axnu, wk2=axlnp)
	aylnp_coef = lnp_test(t0a[agptr], ady, /double, wk1=aynu, wk2=aylnp)
	aflnp_coef = lnp_test(t0a[agptr], (fa[agptr]-mb)/mb, /double, wk1=afnu, wk2=aflnp)	
	caflnp_coef = lnp_test(t0a[agptr], (cfa-ma)/ma, /double, wk1=cafnu, wk2=caflnp)

	axlnp_ptr = where(axlnp eq max(axlnp))
	axnu_max = axnu[axlnp_ptr]
	aylnp_ptr = where(aylnp eq max(aylnp))
	aynu_max = aynu[aylnp_ptr]	
	print, 'Old period = ' + strn(2. / (axnu_max + aynu_max))	
	
	sxb = sxb[bgptr]
	syb = syb[bgptr]
;	scargle,t0b[bgptr],sxb[bgptr],om,bxpx,nu=bnu

; Plot flux periodograms
	pxt = 'Frequency (1/hour)'
	pyt = 'Flux Power (arbitrary units)'
	ptt = 'Periodogram of fluxes for ' 
	dohard = 1
	file = 'fp_r40838144'
	setup_plot, dohard, file
	plot, bfnu, bflnp, PSYM=-3, XTITLE=pxt, YTITLE=pyt, TITLE=ptt + 'r40838144', /XLOG, /YLOG
	oplot, cbfnu, cbflnp, PSYM=-3, COLOR=2
	legend, ['Post-change Raw', '2nd order fit'], PSYM=[1, 1], COLOR=[0, 2], /RIGHT
;	goo = tvrd(true=1)
;	q = color_quan(goo, 1, r, g, b)
;	write_png, 'fp_r40838144.png', goo, r, g, b
	close_plot, dohard, file

	file = 'fp_r40837888'	
	setup_plot, dohard, file
	plot, afnu, aflnp, PSYM=-3, XTITLE=pxt, YTITLE=pyt, TITLE=ptt + 'r40837888', /XLOG, /YLOG
	oplot, cafnu, caflnp, PSYM=-3, COLOR=2
	legend, ['Pre-change Raw', '2nd order fit'], PSYM=[1, 1], COLOR=[0, 2], /RIGHT
	close_plot, dohard, file	
;	goo = tvrd(true=1)
;	q = color_quan(goo, 1, r, g, b)
;	write_png, 'fp_r40837888.png', goo, r, g, b	
	stop

free_lun, lun
; From SNIRAC estimate using f2 = 1.1857685 Jy, b = 0.0064645017 MJy/sr, N/S = 0.000238292
; the ration of the N(snirac)  (8./4196.53  * 1.1857685) / sqrt(bmom[1]) = 1.04 -- reached the 
; Poisson limit, at least I think we did for both experiments in this case
	stop
return
end

pro avg_stack, vec, db, rs, ds, mw, sw
; GOO need to revise for current purposes!!!
	nvec = n_elements(vec)

; Average progressively larger windows, calculate mean and dispersion for each one,
; find all windows such that n is modulo window
	n4 = nvec / 4
; Start with all possible bins
	db = indgen(n4) + 2
	bs = nvec / db
; Find bin sizes that are evenly divisible	
	nrem = nvec mod db
	ptr = where(nrem eq 0, count)
	if (count eq 0) then message, 'No windows are even divisors of sample'
; Keep the evenly divisible bins
	bs = bs[ptr]
	db = db[ptr]
	nwins = n_elements(bs)
	mw = dblarr(nwins)
	sw = dblarr(nwins)
	rs = sw
	ds = sw
	
;	bvec = dblarr(bs[0], nwins) + !VALUES.D_NAN

; Loop over all window sizes, smoothing data and calculating the mean and dispersion.

	for k = 0, nwins-1 do begin
		tvec = dblarr(bs[k])
		for l = 0, bs[k]-1 do tvec[l] = total(vec[l*db[k]:(l+1)*db[k]-1]) / db[k]
;		bvec[0:bs[k]-1, k] = tvec
		mom = moment(tvec)
		mw[k] = mom[0]
		sw[k] = sqrt(mom[1])
	endfor
	rs = sw / sw[0]
	ds = rs - 1. / sqrt(db / db[0])

return
end


pro plot_noise_vs_bin, adb, asw, bdb, bsw
	tek_color
	!P.BACKGROUND=1
	!P.COLOR=0
	!P.PSYM=-3
	!P.THICK=1.25
	!P.CHARSIZE=1.5
	!P.SYMSIZE=1.75
	index = findgen(20000)+2

	window, /free, retain=2, xsize=800, ysize=500

; Plot sample noise vs. binning for one pixel
	xt = 'Bin size (seconds)'
	yt = 'sqrt(variance) [Jy]'
	tt = 'Noise vs. Bin Size'
	plot, adb*0.4, asw, PSYM=1, XTITLE=xt, YTITLE=yt, TITLE=tt, $
	      yrange=[1.E-05, 0.002], xrange=[0.5, 2500.], /YLOG, /XLOG
	oplot, index*0.4, asw[0] / sqrt(index / index[0]), PSYM=-3, COLOR=0
	oplot, bdb*0.4, bsw, PSYM=1, COLOR=3
	oplot, index*0.4, bsw[0] / sqrt(index / index[0]), PSYM=-3, COLOR=3
	!P.PSYM=-3
	legend, ['Pre', 'Post'], PSYM=[1, 1], COLORS=[0, 3], /RIGHT, /TOP, /BOX

	goo = tvrd(true=1)
	q = color_quan(goo, 1, r, g, b)
	write_png, 'noise_coadd.png', goo, r, g, b
	
	set_ps_values, 'noise_coadd.eps', et, cg, mu, sig, xtime, /landscape
	plot, adb*0.4, asw, PSYM=1, XTITLE=xt, YTITLE=yt, TITLE=tt, $
	      yrange=[1.E-05, 0.002], xrange=[0.5, 2500.], /YLOG, /XLOG
	oplot, index*0.4, asw[0] / sqrt(index / index[0]), PSYM=-3, COLOR=0
	oplot, bdb*0.4, bsw, PSYM=1, COLOR=3
	oplot, index*0.4, bsw[0] / sqrt(index / index[0]), PSYM=-3, COLOR=3
	!P.PSYM=-3
	legend, ['Pre', 'Post'], PSYM=[1, 1], COLORS=[0, 3], /RIGHT, /TOP, /BOX
    end_ps
	
	
; Make histograms of % deviation from sqrt(n) noise variation	
;	bsize = 0.001
;	bmax = 0.2
;	bmin = -0.2
;	xt = '% from predict'
;	yt = 'Samples'
;;	tt = 'Noise scaling vs. predict, all pixels, all subframes, AOR 23350272 ch2'
;	for i = 1, n_elements(db)-1 do begin
;		dh = histogram(ds[i, *, *, *], MIN=bmin, MAX=bmax, BINSIZE=bsize, LOCATION=loc)
;		loc = (loc + bsize / 2.0) * 100.
;		if (i eq 1) then plot, loc, dh, XRANGE=[-20.2, 20.2], YRANGE=[0., 2000.], xtitle=xt, $
;		                       ytitle=yt, title=tt, psym=10 $
;		else oplot, loc, dh, psym=10
;		if (i eq 2) then xyouts, -10., 1300., 'DB = ' + strn(db[i]), ALIGNMENT=0.5, ChARSIZE=2.0
;	endfor
;	
;	goo = tvrd(true=1)
;	q = color_quan(goo, 1, r, g, b)
;	write_png, 'noise_vs_model_all.png', goo, r, g, b

return
end

function eclipse_func, x, p
; p is parameters, f0, d, t1, t3, dt
; t is time
    n = n_elements(x)
    model = fltarr(n) + p[0]

    t2 = p[2] + p[4]
    t4 = p[3] + p[4]

; Out of eclipse
;    ptr = where(x le p[2], count)
;    if (count gt 0) then model[ptr] = p[0]

;    ptr = where(x gt t4, count)
;    if (count gt 0) then model[ptr] = p[0]

; Beginning eclipse
    ptr = where(x gt p[2] and x le t2, count)
    if (count gt 0) then $
        model[ptr] = model[ptr] - p[1] * (1.0 - (t2 - x[ptr]) / p[4])

; In eclipse
    ptr = where(x gt t2 and x le p[3], count)
    if (count gt 0) then model[ptr] = model[ptr] - p[1]

; Leaving eclipse
    ptr = where(x gt p[3] and x le t4, count)
    if (count gt 0) then $
        model[ptr] = model[ptr] + p[1] * ((x[ptr] - p[3]) / p[4] - 1.0)

return, model
end

pro fpa1_model2, p, dx, dy, t, model
    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy
    dt = t - min(t)
    if (n_elements(p) ge 7) then scale2 = 1.0 + p[6] * dt else scale2 = 1.0
    if (n_elements(p) eq 8) then scale2 = scale2 + p[7] * dt * dt
    
    model = (1. + scale) * scale2
return
end

function fpa1_xfunc2, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy 
    dt = t - min(t)
    if (n_elements(p) ge 7) then scale2 = 1.0 + p[6] * dt else scale2 = 1.0
    if (n_elements(p) eq 8) then scale2 = scale2 + p[7] * dt * dt
    
    model = model * (1. + scale) * scale2
    model = (flux - model) / err

return, model
end

pro fpa1_model4, p, dx, dy, t, model
    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy + p[6] * dx * dx * dx + p[7] * dx * dx * dy + $
            p[8] * dx * dy * dy + p[9] * dy * dy * dy + $
            p[10] * dx * dx * dx * dx + p[11] * dx * dx * dx * dy + $
            p[12] * dx * dx * dy * dy + p[13] * dx * dy * dy * dy + $
            p[14] * dy * dy * dy * dy
    dt = t - min(t)
    if (n_elements(p) gt 15) then scale2 = 1.0 + p[15] * dt else scale2 = 1.0
    if (n_elements(p) gt 16) then scale2 = scale2 + p[16] * dt * dt
    
    model = (1. + scale) * scale2
return
end

function fpa1_xfunc4, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy + p[6] * dx * dx * dx + p[7] * dx * dx * dy + $
            p[8] * dx * dy * dy + p[9] * dy * dy * dy + $
            p[10] * dx * dx * dx * dx + p[11] * dx * dx * dx * dy + $
            p[12] * dx * dx * dy * dy + p[13] * dx * dy * dy * dy + $
            p[14] * dy * dy * dy * dy
    dt = t - min(t)
    if (n_elements(p) gt 15) then scale2 = 1.0 + p[15] * dt else scale2 = 1.0
    if (n_elements(p) gt 16) then scale2 = scale2 + p[16] * dt * dt
    
    model = model * (1. + scale) * scale2
;    model = (flux-model) / err

return, model
end

function fpa1_func4, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, d, t1, t3, dt, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy + p[6] * dx * dx * dx + p[7] * dx * dx * dy + $
            p[8] * dx * dy * dy + p[9] * dy * dy * dy + $
            p[10] * dx * dx * dx * dx + p[11] * dx * dx * dx * dy + $
            p[12] * dx * dx * dy * dy + p[13] * dx * dy * dy * dy + $
            p[14] * dy * dy * dy * dy
    dt = t - min(t)
    if (n_elements(p) gt 15) then scale2 = 1.0 + p[15] * dt else scale2 = 1.0
    if (n_elements(p) gt 16) then scale2 = scale2 + p[16] * dt * dt
    
    model = model * (1. + scale) * scale2
    model = (flux-model) / err

return, model
end

pro quick_plot2, x, y1, y2, xt, yt, tt, lstr, file
    xmin = min([y1, y2])
    xmax = max([y1, y2])
    xr = xmax - xmin
    xmin = xmin - 0.03 * xr
    xmax = xmax + 0.03 * xr
    plot, x, y1, xtitle=xt, ytitle=yt, title=tt, yrange=[xmin, xmax], PSYM=3
    oplot, x, y2, color=3, PSYM=3
    legend, lstr, colors=[0, 3], PSYM=-3
    goo = tvrd(true=1)
    q = color_quan(goo, 1, r, g, b)
    write_png, file, goo, r, g, b
return
end

pro quick_plot3, x, y1, y2, y3, xt, yt, tt, lstr, file
    xmin = min([y1, y2, y3])
    xmax = max([y1, y2, y3])
    xr = xmax - xmin
    xmin = xmin - 0.03 * xr
    xmax = xmax + 0.03 * xr
    plot, x, y1, xtitle=xt, ytitle=yt, title=tt, yrange=[xmin, xmax], PSYM=3
    oplot, x, y2, color=2, PSYM=3
    oplot, x, y3, color=3, PSYM=3
    legend, lstr, colors=[0, 2, 3], PSYM=-3
    goo = tvrd(true=1)
    q = color_quan(goo, 1, r, g, b)
    write_png, file, goo, r, g, b
    
return
end

pro set_ps_values, fname, et, cg, mu, sig, xtime, LANDSCAPE=landscape
	if (keyword_set(LANDSCAPE)) then begin
		xl = 1
		xp = 0
		xoff = 0.5
		yoff = 10.75
		xsz = 11
		ysz = 7.
	endif else begin
		xl = 0
		xp = 1
		xoff = 0.5
		yoff = 0.5
		xsz = 9.
		ysz = 9.
	endelse
	set_plot, 'PS'
    tek_color
    !P.BACKGROUND=1
    !P.COLOR=0
    !P.PSYM=-3

	!P.THICK=6.
    !P.CHARSIZE=2
    !P.CHARTHICK=1.25
    !P.SYMSIZE=6.
    !X.THICK=1.5
    !Y.THICK=1.5
    et = 5.0
    set_plot,'ps'
    !P.FONT = 0
    cg = 0
    mu = string(181B)
    sig = '!9' + string(115B) + '!6'
    xtime = '!9' + string(180B) + '!6'
    !P.FONT = 0

    device,file=fname,xs=xsz,ys=ysz,/inches,portrait=xp,landscape=xl,xo=xoff,yo=yoff,$
               /helvetica, /isolatin1,/encapsulated, /color

return
end

pro end_ps
	tek_color
	!P.BACKGROUND=1
	!P.COLOR=0
	!P.PSYM=-3
	!P.THICK=1.25
	!P.CHARSIZE=1.5
	!P.SYMSIZE=1.75
	
	device, /CLOSE
	set_plot, 'X'
	!P.FONT=-1   

return
end


pro setup_plot, dohard, file
	if dohard then begin
        tek_color
        !P.BACKGROUND=1
        !P.COLOR=0
        !P.PSYM=-3
        !P.THICK=6.
        !P.CHARSIZE=1.75
        !P.CHARTHICK=1.5
        !P.SYMSIZE=5.
        !X.THICK=1.8
        !Y.THICK=1.8
        set_plot,'ps'
        device,file=file+'.eps',xs=9,ys=9,/inches,/portrait,xo=.5,yo=.5,$
                /TIMES, /isolatin1,/encapsulated, color=1
        !P.FONT = 0
        cg = 0
    endif else begin
    	!P.FONT=-1
		!P.BACKGROUND=1
		!P.COLOR=0
		!P.PSYM=-3
		!P.THICK=1.25
		!P.CHARSIZE=1.5
		!P.SYMSIZE=1.75    
        cg = 0
        window,/FREE,RETAIN=2,xs=500,ys=500
    endelse		
return
end

pro close_plot, dohard, file
	if (dohard) then begin
		device, /CLOSE
		!P.FONT=-1
		set_plot, 'X'
	endif else begin
		goo = tvrd(TRUE=1)
		q = color_quan(goo, 1, r, g, b)
		write_png, file + '.png', goo, r, g, b
	endelse    
return
end