pro shi_m87

ch1_signal = 616. ;mJy
ch1_signal = ch1_signal *1E-3*1E-23
print, 'ch1_signal', ch1_signal, alog10(ch1_signal)
ch1_abmag = flux_to_magab(ch1_signal)
ch1_vegamag = ch1_abmag -2.7 
print, 'ch1 ab, vega', ch1_abmag, ch1_vegamag

ch2_signal = 365. ; mJy
ch2_signal = ch2_signal *1E-3*1E-23
ch2_abmag = flux_to_magab(ch2_signal)
ch2_vegamag = ch2_abmag -3.25
print, 'ch2 ab, vega', ch2_abmag, ch2_vegamag

end
