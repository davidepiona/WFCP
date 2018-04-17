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
replot \
	'plot/plot3.dat' using 1:2 with lines lc rgb "#FFFF00" lw 2 title "Cable 4",\
	'plot/plot3.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot3.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot4.dat' using 1:2 with lines lc rgb "#00FFFF" lw 2 title "Cable 5",\
	'plot/plot4.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot4.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot5.dat' using 1:2 with lines lc rgb "#FF00FF" lw 2 title "Cable 6",\
	'plot/plot5.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot5.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot6.dat' using 1:2 with lines lc rgb "#C0C0C0" lw 2 title "Cable 7",\
	'plot/plot6.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot6.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot7.dat' using 1:2 with lines lc rgb "#800000" lw 2 title "Cable 8",\
	'plot/plot7.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot7.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot8.dat' using 1:2 with lines lc rgb "#808000" lw 2 title "Cable 9",\
	'plot/plot8.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot8.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot9.dat' using 1:2 with lines lc rgb "#008000" lw 2 title "Cable 10",\
	'plot/plot9.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot9.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot10.dat' using 1:2 with lines lc rgb "#800080" lw 2 title "Cable 11",\
	'plot/plot10.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot10.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot11.dat' using 1:2 with lines lc rgb "#008080" lw 2 title "Cable 12",\
	'plot/plot11.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot11.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle