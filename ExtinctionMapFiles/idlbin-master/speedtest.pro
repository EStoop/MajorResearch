pro speedtest
ns = 64
  t = systime(/seconds)
  a = fltarr(ns)
  b = fltarr(ns)
  print, systime(/seconds) - t

  t = systime(/seconds)
  a = fltarr(ns)
  b = a
  print, systime(/seconds) - t
  
  t = systime(/seconds)
  x3 = fltarr(ns)* !VALUES.D_NAN
  y3 = fltarr(ns)* !VALUES.D_NAN
  print, systime(/seconds) - t
  
  t = systime(/seconds)
  x3 = replicate(!VALUES.D_NAN, ns)
  y3 = replicate(!VALUES.D_NAN, ns)
  print, systime(/seconds) - t
  
  t = systime(/seconds)
  x3 = replicate(!VALUES.D_NAN, ns)
  y3 = x3
  print, systime(/seconds) - t

END
