function nearest_neighbors_DT,x,y,chname, DISTANCES=nearest_d,NUMBER=k

  if chname eq 1 then  b = 0.8 else b = 1.0 ; from knutson et al. 2012
  ;help, x
 ; help, y
  ;; Compute Delaunay triangulation
  triangulate,x,y,CONNECTIVITY=c
  help, c
  print, 'n c', n_elements(c), c(n_elements(c) - 25)
  ;limit = n_elements(c)
;  print, 'x', x[0:100]
;  print, 'y', y[0:100]

  n=n_elements(x) 
  nearest=lonarr(k,n,/NOZERO)
  nearest_d=fltarr(k,n,/NOZERO)

;  point = 117980
;  print, 'testing c', c[point], c[point + 1], c[c[point]:c[point+1]-1]
  for point=0L,n-1 do begin 
     ;if point eq 0 then help, point
     ;if point gt limit - 100 then print, 'point', point, c[point], c[point+1]
     ;if point gt 94300 then print, 'point gt 94300', point, c[point], c[point+1] - 1, c[c[point]], c[c[point+1]]
;     nearest[*,point] = !VALUES.F_NAN
;     nearest_d[*,point] = !VALUES.F_NAN

     ;triangulate gets bad results in certain cases ????
     if c[point] gt  c[point+1] - 1 then begin
        ;the problem case
                                ;make up a d that lets me know I should ignore this point
                                ; or just keep the last d which is
                                ; what doing nothing will do I think
                                ; This part seems to happen on the
                                ; transit or eclipse part of the light
                                ; curve, or at least only in those
                                ; AORs where there is a transit or eclipse.
;        print, 'got to the if statement that seems to fail'
 ;       print, point, c[point], c[point+1]
        GOTO, jumpend

      endif  else begin
         p=c[c[point]:c[point+1]-1] ;start with this point's DT neighbors
         d=(x[p]-x[point])^2+((y[p]-y[point])/b)^2
     endelse
                    ; strange triangulate behavior
     for i=1,k do begin 
        s=sort(d)               ; A wasteful priority queue, but anyway
        p=p[s] & d=d[s]
        nearest[i-1,point]=p[0] ; Happy Day: the closest candidate is a NN
        nearest_d[i-1,point]=d[0]
        if i eq k then continue
        
        ;; Add all its neighbors not yet seen
        new=c[c[p[0]]:c[p[0]+1]-1]
        nnew=n_elements(new) 
        already_got=[p,nearest[0:i-1,point],point]
        ngot=n_elements(already_got) 
        wh=where(long(total(rebin(already_got,ngot,nnew,/SAMPLE) ne $
                            rebin(transpose(new),ngot,nnew,/SAMPLE),1)) $
                 eq ngot, cnt)
        if cnt gt 0 then begin 
           new=new[wh]
           p=[p[1:*],new]
           d=[d[1:*],(x[new]-x[point])^2+((y[new]-y[point])/b)^2]
        endif else begin 
           p=p[1:*]
           d=d[1:*]
        endelse 
     endfor
     JUMPEND: 

  endfor 
  if arg_present(nearest_d) then nearest_d=sqrt(nearest_d)

   return, nearest

end 
