set autoscale

set term wxt title '2.811260e+14'
plot \
	'plot/plot_datastd_1.dat' using 1:2 with lines lc rgb "#00FF00" lw 2 title "Cable 2",\
	'plot/plot_datastd_1.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_1.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_2.dat' using 1:2 with lines lc rgb "#0000FF" lw 2 title "Cable 3",\
	'plot/plot_datastd_2.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_2.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_3.dat' using 1:2 with lines lc rgb "#FFFF00" lw 2 title "Cable 4",\
	'plot/plot_datastd_3.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_3.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_4.dat' using 1:2 with lines lc rgb "#00FFFF" lw 2 title "Cable 5",\
	'plot/plot_datastd_4.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_4.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_5.dat' using 1:2 with lines lc rgb "#FF00FF" lw 2 title "Cable 6",\
	'plot/plot_datastd_5.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_5.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_6.dat' using 1:2 with lines lc rgb "#C0C0C0" lw 2 title "Cable 7",\
	'plot/plot_datastd_6.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_6.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_7.dat' using 1:2 with lines lc rgb "#800000" lw 2 title "Cable 8",\
	'plot/plot_datastd_7.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_7.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_9.dat' using 1:2 with lines lc rgb "#008000" lw 2 title "Cable 10",\
	'plot/plot_datastd_9.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_9.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_11.dat' using 1:2 with lines lc rgb "#008080" lw 2 title "Cable 12",\
	'plot/plot_datastd_11.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_11.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_12.dat' using 1:2 with lines lc rgb "#000080" lw 2 title "Cable 13",\
	'plot/plot_datastd_12.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_12.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_datastd_13.dat' using 1:2 with lines lc rgb "#808080" lw 2 title "Cable 14",\
	'plot/plot_datastd_13.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_datastd_13.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle