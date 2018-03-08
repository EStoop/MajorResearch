pro write_irafspec,spec,$     ; the spectrum vector
                   lambda,$   ; the associated wavelength vector
                   outfile    ; output fits file name as string


length=n_elements(spec)
lambda_start=lambda[0]
lambda_stop=lambda[length-1]
dispersion=(lambda_stop-lambda_start)/length

lambda_linear=findgen(length)
lambda_linear=lambda_start+(lambda_linear*dispersion)


; do the linear interpolation
linterp,lambda,spec,lambda_linear,out




mkhdr,hd,spec ; make the header
sxaddpar,hd,'CTYPE1','WAVE-WAV'
sxaddpar,hd,'CUNIT1','Angstrom'
sxaddpar,hd,'CRVAL1',lambda_start ; value at crpix1
sxaddpar,hd,'CD1_1',dispersion ; linear dispersion
sxaddpar,hd,'CRPIX1',0  ; the first pixel

writefits,outfile,out,hd


end
