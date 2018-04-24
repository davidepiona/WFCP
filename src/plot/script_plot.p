set autoscale
plot \
	'plot/plot0.dat' using 1:2 with lines lc rgb "#FF0000" lw 2 title "Cable 1",\
	'plot/plot0.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot0.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot1.dat' using 1:2 with lines lc rgb "#00FF00" lw 2 title "Cable 2",\
	'plot/plot1.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot1.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot2.dat' using 1:2 with lines lc rgb "#0000FF" lw 2 title "Cable 3",\
	'plot/plot2.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot2.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle