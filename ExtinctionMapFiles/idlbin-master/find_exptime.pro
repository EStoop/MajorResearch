Function find_exptime, source_mJy, ch
  ;;this function takes the flux density
  ;;of a source in mJy as well as the
  ;;channel and makes an estimate of which exposure time is
  ;;appropriate for that source.  The goal was to stay below 2/3 full
  ;;well and above 10-15% full well.

  ;;based on warm mission saturation limits on IRAC website
  ;;http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/warmimgcharacteristics/
  if ch eq '1' then begin
    case 1 of
       (source_mJy ge 2000): exptime = 0.02
       (source_mJy lt 2000) and (source_mJy ge 460): exptime = 0.1
       (source_mJy lt 460) and (source_mJy ge 100): exptime = 0.4
       (source_mJy lt 100) and (source_mJy ge 20): exptime = 2.
       (source_mJy lt 33) and (source_mJy ge 16): exptime = 6.
       (source_mJy lt 16) and (source_mJy ge 7): exptime = 12.
       (source_mJy lt 7) and (source_mJy ge 1.7): exptime = 30.
       (source_mJy lt 1.7) and (source_mJy ge 0.5): exptime = 100.
       (source_mJy lt 0.5): begin
          print, "Not a candidate for Spitzer staring observations"
          return, 0
       end
    endcase
 endif

  if ch eq '2' then begin
    case 1 of
       (source_mJy ge 2400): exptime = 0.02
       (source_mJy lt 2400) and (source_mJy ge 600): exptime = 0.1
       (source_mJy lt 600) and (source_mJy ge 120): exptime = 0.4
       (source_mJy lt 120) and (source_mJy ge 22): exptime = 2.
       (source_mJy lt 40) and (source_mJy ge 20): exptime = 6.
       (source_mJy lt 20) and (source_mJy ge 8): exptime = 12.
       (source_mJy lt 8) and (source_mJy ge 2): exptime = 30.
       (source_mJy lt 2) and (source_mJy ge 0.5): exptime = 100.
       (source_mJy lt 0.5): begin
          print, "Not a candidate for Spitzer staring observations"
          return, 0
       end
    endcase
 endif

return, exptime
end


