function splitsimbadstring, firstline, starname = starname, ra = ra, dec = dec, pmra = pmra, pmdec = pmdec
     ;;firstline = '1|      0.99|KOI-157.02               |Pl |297.1150-41.9092           |            ~ ~ |     ~|     ~|     ~|     ~|     ~|~              |  67|   1'
  firstsplit = firstline.Split('\|')
     starname = firstsplit[2]
     starname = starname.Compress()
     pos = firstsplit[4]
     pos = pos.Compress()
     pm = firstsplit[5]         ;mas/yr
     pm = pm.Trim()
     minus = pos.Contains('-')
     if minus gt 0 then begin
        negdec = pos.Split('\-')
        ra = float(negdec[0])
        dec = float(negdec[1])
        dec = dec*(-1)
     endif else begin
        posdec = pos.Split('\+')
        ra = float(posdec[0])
        dec = float(posdec[1])
     endelse
     if pm ne '~ ~' then begin
        both = pm.Split(' ')
        pmra = both[0]
        pmdec = both[1]
     endif else begin  ;;unknown proper motions
        pmra = 0
        pmdec = 0
     endelse

     return, 0
  end

  
