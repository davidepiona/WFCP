set autoscale
plot \
	'plot/plot_datastd_0.dat' using 1:2 with lines lc rgb "#FF0000" lw 2 title "Cable 1",\
	'plot/plot_datastd_0.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_0.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_1.dat' using 1:2 with lines lc rgb "#00FF00" lw 2 title "Cable 2",\
	'plot/plot_datastd_1.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_1.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_2.dat' using 1:2 with lines lc rgb "#0000FF" lw 2 title "Cable 3",\
	'plot/plot_datastd_2.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_2.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle