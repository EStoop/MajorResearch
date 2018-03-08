pro box_centroider, image, sigma2, xmax, ymax, boxwidth, backboxwidth, boxborder, $
                    x0, y0, f0, b, xs, ys, fs, bs, c, cb, np, xfwhm, yfwhm, MMM=mmm
;+
; NAME:
;    BOX_CENTROIDER
; 
; PURPOSE:
;    Calculates centroids for a source using a simple 1st moment box centroider
; 
; INPUTS:
;    IMAGE: 2d float/double array containing input image, input image is presumed to 
;           be have bad pixels NaNed out
;    SIGMA2: 2d float/double array containing estimate of square of per pixel image 
;            uncertainty
;    XMAX: float/double scalar containing x position of peak pixel to centroid around
;    YMAX: float/double scalar containing y position of peak pixel to centroid around
;    BOXWIDTH: integer scalar, for size of centroiding box, actually halfboxwidth - 0.5
;    BACKBOXWIDTH: integer scalar, for size of rectangular background box
;    BOXBORDER: integer scalar, offset of background annulus in pixels from box aperture
;
; CONTROL KEYWORD:
;    MMM: if set, then backround estimate will be calculated using mmm.pro
;
; OUTPUTS:
;    X0: double scalar X position of centroid
;    Y0: double scalar Y position of centroid
;    F0: double scalar flux of source as summed in box after background removal
;    B: double scalar background level
;    XS: double scalar uncertainty in X centroid position
;    YS: double scalar uncertainty in Y centroid position
;    FS: double scalar uncertainty in measured flux accounting for uncertainty in background
;    BS: double scalar uncertainty in background
;    C: fixed scalar number of good pixels in aperture
;    CB: fixed scalar number of good pixels in background aperture
;    NP: fixed scalar number of noise pixels
;    XFWHM: full width at half maximum in x-coordinate of centroid
;    YFWHM: full width at half maximum of y-coordinate of centroid
;
; ALGORITHM:
;    
; HISTORY:
;    Added source width calculation 27 July 2011 SJC
;    Recoding from get_centroids_nf procedure 26 May 2010 SJC
;    
;-

; Set sigma clipping threshold
  sigma = 3

; Get image size
  sz = size(image)
  nx = sz[1]
  ny = sz[2]
    
; Arrays holding x and y coordinates of each pixel -- allow arbitrary sized 
; images
  xx = findgen(nx) # replicate(1.0, ny)
  yy = replicate(1.0, nx) # findgen(ny)
  
; Set box limits
  xa = xmax-boxwidth > 0 < (nx-1)
  xb = xmax+boxwidth < (nx-1) > 0
  ya = ymax-boxwidth > 0 < (ny-1)
  yb = ymax+boxwidth < (ny-1) > 0

; Set background regions
  xamin = xa - backboxwidth - boxborder > 0 < xa
  xamax = xa - boxborder > 0 < xa
  xbmax = xb + backboxwidth + boxborder > xb < (nx-1)
  xbmin = xb + boxborder > xb < (nx-1)
  yamin = ya - backboxwidth - boxborder > 0 < ya
  yamax = ya - boxborder > 0 < ya
  ybmax = yb + backboxwidth + boxborder > yb < (ny-1)
  ybmin = yb + boxborder > yb < (ny-1)

; Segment image into box and everything else
  mask = image * 0
  mask[xa:xb, ya:yb] = 1
  bptr = where(xx le xamax and xx ge xamin, bcount)
  if (bcount gt 0) then mask[bptr] = -1
  bptr = where(xx le xbmax and xx ge xbmin, bcount)
  if (bcount gt 0) then mask[bptr] = -1
  bptr = where(yy le yamax and yy ge yamin, bcount)
  if (bcount gt 0) then mask[bptr] = -1
  bptr = where(yy le ybmax and yy ge ybmin, bcount)
  if (bcount gt 0) then mask[bptr] = -1

; remove background - median of non box region
  bptr = where(mask eq -1 and image eq image, bcount)
  if (bcount gt 0) then begin 
    if (keyword_set(MMM)) then begin
      mmm, image[bptr], back, bsig
      bsig = bsig * bsig
    endif else begin
      back = median(image[bptr]) 
; One pass of outlier rejection for the background, could use mmm as well
      bmom = moment(image[bptr])
      bsig = bmom[1]
      bbptr = where(abs(image[bptr] - back) le sigma * sqrt(bsig), $
                         bbcount)
      if (bbcount gt 0) then begin
        back = median(image[bptr[bbptr]])
        bmom = moment(image[bptr[bbptr]] - back)
        bsig = bmom[1]
        bcount = bbcount
      endif else begin
        back = 0.0D
        bsig = 0.0D
      endelse
    endelse
  endif else begin
    back = 0.0D
    bsig = 0.0D
  endelse

; Perform background subtraction
  image = image - back
  b = back
  
  npixels_in_box = (2*boxwidth+1) * (2*boxwidth+1)

; Perform centroid in 2*boxwidth+1 pixel box centered on peak pixel
  gptr = where(mask eq 1 and image eq image, gcount)
  if (gcount gt 0.5 * npixels_in_box) then begin
    fluxsum = total(image[gptr])
    flux2sum = total(image[gptr] * image[gptr])
    fluxsum2 = fluxsum * fluxsum
    x0 = total(xx[gptr] * float(image[gptr])) / fluxsum
    y0 = total(yy[gptr] * float(image[gptr])) / fluxsum
    f0 = fluxsum
; Now calculate the variance 
; note for normal distributions, fwhm = 2*sqrt(2 * ln 2) sigma ~ 2.35482 sigma
    dx = xx[gptr] - x0
    dy = yy[gptr] - y0
    xfwhm = total(dx * dx * float(image[gptr])) / fluxsum
    yfwhm = total(dy * dy * float(image[gptr])) / fluxsum
    sigscale = 2.D * sqrt(2.D * alog(2.D))
    xfwhm = sigscale * sqrt(xfwhm)
    yfwhm = sigscale * sqrt(yfwhm)
    
; Use summed flux as measure of total flux and calculate noise pixels   
    np = flux2sum / fluxsum
; Number of pixels in source and background regions
    c = float(gcount)
    cb = float(bcount)
    
; If a bad pixel exists in the source aperture return a flux of NaN.
    xptr = where(mask eq 1 and image ne image, xcount)
    if (xcount gt 0) then f = !VALUES.D_NAN else f = fluxsum
  
; Uncertainties
    xs = total(dx * dx * sigma2[gptr])
    ys = total(dy * dy * sigma2[gptr])
    xs = xs / fluxsum2
    ys = ys / fluxsum2
    bs = c * c * bsig / cb
    fs = total(sigma2[gptr]) + bs
  endif else begin
    print, 'BOX_CENTROIDER: Less than 50% good pixels in source box'
    x0 = !VALUES.D_NAN
    y0 = !VALUES.D_NAN
    f0 = !VALUES.D_NAN
    b = !VALUES.D_NAN
    xs = !VALUES.D_NAN
    ys = !VALUES.D_NAN
    bs = !VALUES.D_NAN
    fs = !VALUES.D_NAN
    xfwhm = !VALUES.D_NAN
    yfwhm = !VALUES.D_NAN   
  endelse

return
end

