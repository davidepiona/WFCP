set autoscale

<<<<<<< HEAD
set term wxt title '1.000000e+20'
=======
<<<<<<< HEAD
set term wxt title '1.000000e+20'
=======
set term wxt title '2.000166e+11'
plot \
	'plot/plot_datastd_0.dat' using 1:2 with lines lc rgb "#FF0000" lw 2 title "Cable 1",\
	'plot/plot_datastd_0.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_0.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_1.dat' using 1:2 with lines lc rgb "#00FF00" lw 2 title "Cable 2",\
	'plot/plot_datastd_1.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_1.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
>>>>>>> b1dde4eb99f6b59940d36baaaf9add85c318a20f
>>>>>>> a07193c702de53b440839013afe899e0139e5c77
