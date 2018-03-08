pro col_define_example
  
  twopi=6.2831853071795865D
  i = indgen(256)
  x = cos(float(i)/256*twopi) &  y=sin(float(i)/256*twopi)

  !p.multi = [0,2,2]
  
  col_define,39
  plot,[-1,1],[-1,1],/nodata,col=!col.wheat
  for i=0,255 do oplot,[0,x[i]],[0,y[i]],col=!col.ct[i]

  col_define,39,16
  plot,[-1,1],[-1,1],/nodata,col=!col.cyan
  for i=0,255 do oplot,[0,x[i]],[0,y[i]],col=!col.ct[i]

  col_define,24
  plot,[-1,1],[-1,1],/nodata,col=!col.chartreuse
  for i=0,255 do oplot,[0,x[i]],[0,y[i]],col=!col.ct[i]

  col_define,24,48
  plot,[-1,1],[-1,1],/nodata,col=!col.orange
  for i=0,255 do oplot,[0,x[i]],[0,y[i]],col=!col.ct[i]
  
end

  
