function calc_earthpoint_pitch, bmjd0, pa_earth

  ;;first read in the csv file with the dates and downlink pitch
  ;;angles from Effertz
  ;;pitch_angles  = read_csv_jim('/Users/jkrick/irac_warm/Spitzer_pitch.csv')
  jd = julday(pa_earth.month, pa_earth.day, pa_earth.year)
  mjd = jd - 2400000.5

  
  ;;which value in csv file has the same mjd?
  bin = Value_Locate(mjd, bmjd0) ;;choose the closest one
  pe = pa_earth.pitch_earth(bin)
        
  return, pe  
  

end
