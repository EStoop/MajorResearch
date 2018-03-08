pro make_psfigure

fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits', iracdata, iracheader

ps_open, filename='/Users/jkrick/noao/gemini/clusteridl.ch1.ps',/portrait

plotimage, xrange=[2159,2332],yrange=[1880,2027], bytscl(iracdata, min = 00.04 ,max = 0.5) ,$;
 /preserve_aspect, /noaxes, ncolors=8

ps_close, /noprint,/noid



ps_open, filename='/Users/jkrick/noao/gemini/clusteridl.acs.ps',/portrait

      acsdata = mrdfits('/Users/jkrick/hst/raw/wholeacs.fits', range=[12874, 14619])
      plotimage, xrange=[19381,21465],$
                 bytscl(acsdata, min = -0.01, max = 0.1),$
                 /preserve_aspect, /noaxes, ncolors=8


ps_close, /noprint,/noid

end
