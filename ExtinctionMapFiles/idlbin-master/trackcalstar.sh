#!/bin/tcsh
 #setenv IDL_PATH +/Users/jkrick/idlbin:+/Users/jkrick/idlbin/cm:+/Volumes/Nutella/Users/jkrick/idlbin:+/Volumes/Nutella/opt/idlutils/idlutils-v5_2_0/pro:+/Users/jkrick/idlbin/coyote:+/Users/jkrick/idlbin/idlastro/pro

osascript -e 'tell application "Finder" to open location "nfs://ssceng1/ssc/Http/Docs/Instruments/IRAC/"'
/Applications/exelis/idl85/bin/idl -e  "track_calstar_photometry, '/Users/jkrick/irac_warm/calstars/track_calstar_ch2.sav', 2, /make_plot,/binning"
/Applications/exelis/idl85/bin/idl -e  "track_calstar_photometry, '/Users/jkrick/irac_warm/calstars/track_calstar_ch1.sav', 1, /make_plot,/binning"
