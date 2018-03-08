pro cluster_ellipticals

restore, '/Users/jkrick/idlbin/objectnew.sav'

ellip = [1250,974,1094,1120,1766,2009,2540]

plothist, objectnew[ellip].acsmag - objectnew[ellip].irac1mag, bin = 0.1, xtitle = 'ACS F814W - [3.6]', xrange=[1,5], charthick = 1


end
