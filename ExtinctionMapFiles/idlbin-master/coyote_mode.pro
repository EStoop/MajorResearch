function coyote_mode, array
 array = array[Sort(array)] 
   wh = where(array ne Shift(array,-1), cnt)
   if cnt eq 0 then mode = array[0] else begin
      void = Max(wh-[-1,wh], mxpos)
      mode = array[wh[mxpos]]
   endelse 
   Print, 'coyote', mode   
return, mode
end
