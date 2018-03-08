        PRO multipole, intemps, tcb, sigt, di_amp, sigdp, di_glon, di_glat, $
                       sigdpd, galexc=galexc, residmap=residmap, $
                       badval=badval, weights=weights, quadmap=quadmap, $
                       qparms=qparms,qsigma=qsigma
;+NAME/ONE-LINE DESCRIPTION:
;    MULTIPOLE subtracts dipole + quadrupole terms from a scalar sky map.
;
; DESCRIPTION:
;    An input intensity map (with optional weights) is the basis of
;    a least-squares fit for the position, monopole value, dipole
;    amplitude, quadrupole amplitude and quadrupole coefficients of the 
;    CMBR.  The routine returns these values as well 
;    as a residual map with the dipole plus quadrupole subtracted.  
;    The user can also specify a range of galactic latitude to exclude 
;    from the fit. The residual is subtracted from the whole map. A
;    pure quadrupole map is also created. The routine assumes that the 
;    chi-square per degree of freedom of the fit is 1 and computes errors 
;    on the dipole amplitude and direction.
;
; CALLING SEQUENCE:
;  MULTIPOLE, inmap, monopole, monopole_sig diampl, diampl_sig, diglon, diglat,$
;           digl_sig, [galexc=...], [residmap=residmap], [weights=weights],$
;           [badval=badval],[qparms=qparms],[qsigma=qsigma],[quadmap=quadmap]
;
; ARGUMENTS (I=input, O=output, []=optional):
;    inmap        I    flt arr   2-D input skymap, in right-T format
;    monopole     O    dble      Mean CMBR intensity 
;    monopole_sig O    dble      Sigma of CMBR intensity
;    diampl       O    dlbe      Dipole amplitude
;    diampl_sig   O    dlbe      Sigma of dipole amplitude
;    diglon       O    dlbe      Galactic longitude of dipole pole (deg)
;    diglat       O    dlbe      Galactic latitude of dipole pole (deg)
;    digl_sig     O    dlbe      Sigma of dipole direction (deg)
;   [galexc]     [I]   flt       Absolute value of minimum galactic
;                                      latitude to include in fit (deg),
;                                      default=0.
;   [residmap]   [O]   flt arr   2-D output residual map
;   [weights]    [I]   flt arr   2-D map of weights corresponding to
;                                       input map, default=uniform wts
;   [badval]     [I]   flt       Flag value for bad pixels, default=0.
;   [qparms]     [O]   dlbe      Quadrupole coefficients
;   [qsigma]     [O]   dlbe      Sigma's of quadrupole coefficients
;   [quadmap]    [O]   dlbe      2-D output of quadrupole map
;
; WARNINGS: 
;    1.  Only a full sky right-T can be used as input.  Pixelization
;           is assumed to be in ecliptic coordinates!
;    2.  It is generally wise to specify a galactic latitude exclusion
;           of 20 degrees or more for believable results.
;    3.  Flagged pixels will show up as zeros in the residual map.
;
; EXAMPLE:
;    If "inmap" is a map of temperatures in which empty pixels have
;    a value of zero then a simple call for an unweighted fit would be
;
;    UIDL> multipole,inmap,avtemp,tsig,diamp,ampsig,lon,lat,dir_sig,$
;           galexc=20,qparms=qparms
;    UIDL> print,avtemp,tsig,diamp,ampsig,lon,lat,dir_sig,qparms
;             2.7256  0.000005  0.003345  0.000016  265.60  48.34  0.24
;             -0.0368 0.0254 -0.0266 -0.0196 0.0120 0.0278
;   If, e.g., "residmap=residual" had also been included in the call 
;   then the variable "residual" would contain a residual skymap after
;   subtracting a 2.7256 K average plus a 3.345 mK dipole plus a 
;   quadrupole from inmap.
;
;   If in addition to inmap we have a second map, e.g. "sigmas" that
;   contains the weight at each pixel, then we could also have
;   included ", weights=sigmas" in the call.
;#
; COMMON BLOCKS:  none
;
; PROCEDURE:
;    There is a fair amount of coordinate conversion done, so the first
;    thing that is done is to determine the resolution of the skymap
;    from its dimensions.  The procedure is:
;       (1)  Determine map resolution
;       (2)  "De-rasterize" the data into a 1-D pixel list
;       (3)  Cull pixels with flag values from the list
;       (4)  Cull low galactic latitude pixels from the list
;       (5)  Convert pixel numbers of remaining pixels to unit vectors
;       (6)  Using unit vectors, do a SVD fit to the dipole +quadrupole in
;               Cartesian coordinates, i.e. Model = T0 + Td.u + Tq.uq, where
;               the second term is the dot product of the dipole (of
;               amplitude Td) and the unit vector u and the third term
;               is the dot product of the quadrupole (of amplitude Tq) and
;               the Legendre P2 unit vector u.
;       (7)  Compute sigma's using output from SVD matrices.
;       (8)  Convert the best-fit unit vector into galactic coords (for dipole)
;       (9)  Compute the sigma for the direction by recomputing the unit vector
;               for the dipole + sigmas and for the dipole - sigmas.
;      (10)  Create a dipole + quadrupole map and subtract from the input map
;               if a residual map is desired
;      (11)  Create a quadrupole model map if desired.
;
; LIBRARY CALLS:
;     coorconv
;     pix2xy
;     xy2pix
;
; MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  May 1992
;    Modified by Celine Groden, USRA,  December 1992.  SPR # 10327
;    Modified by Gene Eplee, General Sciences Corp, Dec 1992
;        to compute sigmas for fit.
;    SPR 10941  May 14 93  Extensively modified to calculate quadrupole term.
;
;.title
;Routine MULTIPOLE
;-
;  Get cube resolution and verify that map is a scalar
;
        ON_ERROR, 2
        dims = SIZE (intemps)
        cuberes = FIX(3.33 * alog10(dims(1)/4) + 1)
        IF (dims(0) NE 2) THEN MESSAGE, "Map is not 2-D"
;
;  Set defaults for optional parameters (e.g. uniform weighting) 
;
        IF (N_ELEMENTS(badval) EQ 0) THEN badval = 0.
        IF (N_ELEMENTS(galexc) EQ 0) THEN galexc = 0.                       
        IF (N_ELEMENTS(weights) EQ 0) THEN weights = 0.*intemps + 1.
;        
;  Get the 1-D index of all non-flag data points, and create 
;  corresponding 1-D pixel-ordered lists of good data, weights,
;  and the pixel numbers.
;
        gooddata = WHERE (intemps NE badval)
        cbtemps = intemps(gooddata)                        ; good data
        wtvec = weights(gooddata)                          ; weights
        xcoords = gooddata MOD dims(1)
        ycoords = FIX (gooddata/dims(1))
        pixlist = xy2pix (xcoords, ycoords, res=cuberes)   ; pixel nums
;
;  Get the galactic coords of all pixels in the list, and filter
;  out the ones whose latitude is too low.
;
        resstr = 'R' + STRTRIM(STRING(cuberes),2)
        galco = coorconv (pixlist, infmt='p', outfmt='l', inco=resstr, $
                             outco='g')
        highlat = WHERE (ABS(galco(*,1)) GE galexc)
        bsz=SIZE(highlat)
        IF (bsz(0) EQ 0) THEN BEGIN
           PRINT,'There is no good data above specifed galactic latitude'
           RETURN
        ENDIF
        pixlist = pixlist(highlat)                        ; culled pixel list
        cbtemps = cbtemps(highlat)                        ; culled data
        wtvec = wtvec(highlat)                            ; culled weights
        avtemp = total(cbtemps)/n_elements(cbtemps)
        cbtemps = cbtemps - avtemp                        ; create a zero mean
;
;  Create the unit vector for each pixel, then fill the normal matrix
;  QMAT for the fit. The unit vectors are the derivatives of the model, and
;  the fit equation is  QMAT # DPARMS = RHS
;
        uvecs = coorconv (pixlist, infmt='p', outfmt='u', inco=resstr, $
                             outco='g')
        uvecs = DOUBLE(uvecs)
;
npix=SIZE(cbtemps)
K = DBLARR(9,npix(1)) 
K(0,*) = 1.0
K(1,*) = uvecs(*,0)
K(2,*) = uvecs(*,1)
K(3,*) = uvecs(*,2)
K(4,*) = 0.5*(3.*uvecs(*,2)^2. - 1)
K(5,*) = 2.*uvecs(*,2)*uvecs(*,0)
K(6,*) = 2.*uvecs(*,2)*uvecs(*,1)
K(7,*) = (uvecs(*,0)^2. - uvecs(*,1)^2.)
K(8,*) = 2.*uvecs(*,1)*uvecs(*,0)
    
N1 = TOTAL(cbtemps * wtvec * K(0,*))
N2 = TOTAL(cbtemps * wtvec * K(1,*))
N3 = TOTAL(cbtemps * wtvec * K(2,*))
N4 = TOTAL(cbtemps * wtvec * K(3,*))
N5 = TOTAL(cbtemps * wtvec * K(4,*))
N6 = TOTAL(cbtemps * wtvec * K(5,*))
N7 = TOTAL(cbtemps * wtvec * K(6,*))
N8 = TOTAL(cbtemps * wtvec * K(7,*))
N9 = TOTAL(cbtemps * wtvec * K(8,*))
N  = [N1,N2,N3,N4,N5,N6,N7,N8,N9]

M = DBLARR(9,9)
FOR i = 0,8 DO BEGIN
    FOR j = 0,8 DO BEGIN
        M(i,j) = TOTAL( K(i,*)*K(j,*)*wtvec )
    ENDFOR
ENDFOR

; use SVD to invert matrix M and solve...
; ...A is the inverse of the Single Value Decomposition of M

SVD, M, w, U, V
Small = WHERE ( w LT MAX(w) *1.0e-06, count )
IF (Count ne 0) THEN  w(small) = 0.0
Wp = DBLARR(N_elements(w),N_elements(w))
FOR i =0,(N_elements(w) - 1) DO IF(w(i) ne 0) THEN Wp(i,i) = 1./w(i)
A = V # Wp # TRANSPOSE(U)
dparms = V # Wp # (TRANSPOSE(U) # N)

tcb = dparms(0) +avtemp
Q1       = dparms(4)
Q2       = dparms(5)
Q3       = dparms(6)
Q4       = dparms(7)
Q5       = dparms(8)
qrms = SQRT( (4./15.)*(0.75*Q1^2. + Q2^2. + Q3^2. + Q4^2. + Q5^2.) )

; The errors in the fitted coefficients are given by the diag elements of A if
; weighting is used!!

variances  = DBLARR(9)
sigmas  = DBLARR(9)
FOR j = 0,8 DO BEGIN
    variances(j)   = A(j,j)
    sigmas(j)   = SQRT( variances(j) )
ENDFOR
Errq = SQRT( (variances(4)*0.75^2.*Q1^2. + variances(5)*Q2^2. + $
       variances(6)*Q3^2. + $
       variances(7)*Q4^2. + variances(8)*Q5^2.) ) * (4./(15.*qrms))

;  Normalize the cartesian representation of the dipole to get the unit
;  vector of its direction, and convert to galactic coordinates.  The 
;  squared modulus of the three componenets is the dipole amplitude.
;
        di_amp = SQRT (TOTAL (dparms(1:3)^2))
        di_uvec = dparms(1:3)/di_amp
        galcoords = coorconv(di_uvec, infmt='u', outfmt='l')
        di_glon = galcoords(0)
        di_glat = galcoords(1)
; Parameters and errors for Quadrupole
qparms=[Q1,Q2,Q3,Q4,Q5,qrms]
qsigma=[sigmas(4),sigmas(5),sigmas(6),sigmas(7),sigmas(8),Errq]
;
;  Error in monopole amplitude.
;
	sigt = sigmas(0)
;
;  Error in dipole amplitude.
;
	sigdp = sqrt(variances(1)+variances(2)+variances(3))
;
;  Compute the error in dipole direction by recomputing the direciton for
;      the original unit vector plus sigma and for the original unit vector
;      minus sigma.
;
	sigdpu = dblarr(3)
	sigdpu(0) = dparms(1) + sigmas(1)
	sigdpu(1) = dparms(2) + sigmas(2)
	sigdpu(2) = dparms(3) + sigmas(3)
        sigd_amp = SQRT(TOTAL(sigdpu^2))
        sigd_uvec = sigdpu/sigd_amp
        sigcoord = coorconv(sigd_uvec, infmt='u', outfmt='l')
        sigh_glon = abs(sigcoord(0) - di_glon)
        sigh_glat = abs(sigcoord(1) - di_glat)

	sigdpu(0) = dparms(1) - sigmas(1)
	sigdpu(1) = dparms(2) - sigmas(2)
	sigdpu(2) = dparms(3) - sigmas(3)
        sigd_amp = SQRT(TOTAL(sigdpu^2))
        sigd_uvec = sigdpu/sigd_amp
        sigcoord = coorconv(sigd_uvec, infmt='u', outfmt='l')
        sigl_glon = abs(sigcoord(0) - di_glon)
        sigl_glat = abs(sigcoord(1) - di_glat)
;
;  The average difference of the dipole directions is error in the direction.
;
	sigdpd = (sigh_glon + sigh_glat + sigl_glon + sigl_glat)/4.0d

;
;  We now have all fit parameters, so create an all-sky map holding only
;  the dipole.  First generate a complete pixel list holding pixel
;  numbers 0-(size of cube), and create the unit vector for each point.
;
        npix = 3 * 2^(2*cuberes-1) 
        allpix = lindgen(npix)
        alluvecs = coorconv (allpix, infmt='p', outfmt='u', $
                                     inco=resstr, outco='g')
        QK=DBLARR(5,npix)
        QK(0,*) = 0.5*(3.*alluvecs(*,2)^2. - 1)
        QK(1,*) = 2.*alluvecs(*,2)*alluvecs(*,0)
        QK(2,*) = 2.*alluvecs(*,2)*alluvecs(*,1)
        QK(3,*) = (alluvecs(*,0)^2. - alluvecs(*,1)^2.)
        QK(4,*) = 2.*alluvecs(*,1)*alluvecs(*,0)
;
;  Now create a model sky by taking the dot product of each pixel's
;  unit vector with the dipole, and adding in the monopole.
;
        tsky_model = fltarr(npix)
        for j=0,npix-1 do tsky_model(j) = tcb + $
                    total(alluvecs(j,*) * dparms(1:3)) 
        di_model=tsky_model
        for j=0,npix-1 do tsky_model(j)=tsky_model(j) + $
                    total(qk(*,j)*dparms(4:8))
;                    (Q1*QK(0,*)) + (Q2*QK(1,*)) + (Q3*QK(2,*)) + $
;                    (Q4*QK(3,*)) + (Q5*QK(4,*))

;  Rasterize the model sky, and flag as bad any pixels that are 
;  bad in the original data. Then difference it with the original data
;  to create a residual map.  (Residuals will be zero at flagged points.)
;
        pix2xy, allpix, res=cuberes, data=tsky_model, raster=outmap
        outmap(where(intemps eq badval)) = badval
        residmap = intemps - outmap
        pix2xy, allpix, res=cuberes, data=di_model, raster=dimap
        dimap(where(intemps eq badval)) = badval
        quadmap = outmap-dimap
        end

;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.


