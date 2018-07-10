set autoscale

set term wxt title '9.154860e+06'
plot \
	'plot/plot_data_17_1.dat' using 1:2 with lines lc rgb "#00FF00" lw 2 title "Cable 2",\
	'plot/plot_data_17_1.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_1.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_2.dat' using 1:2 with lines lc rgb "#0000FF" lw 2 title "Cable 3",\
	'plot/plot_data_17_2.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_2.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_3.dat' using 1:2 with lines lc rgb "#FFFF00" lw 2 title "Cable 4",\
	'plot/plot_data_17_3.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_3.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_4.dat' using 1:2 with lines lc rgb "#00FFFF" lw 2 title "Cable 5",\
	'plot/plot_data_17_4.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_4.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_5.dat' using 1:2 with lines lc rgb "#FF00FF" lw 2 title "Cable 6",\
	'plot/plot_data_17_5.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_5.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_6.dat' using 1:2 with lines lc rgb "#C0C0C0" lw 2 title "Cable 7",\
	'plot/plot_data_17_6.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_6.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_7.dat' using 1:2 with lines lc rgb "#800000" lw 2 title "Cable 8",\
	'plot/plot_data_17_7.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_7.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
replot \
	'plot/plot_data_17_8.dat' using 1:2 with lines lc rgb "#808000" lw 2 title "Cable 9",\
	'plot/plot_data_17_8.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
	'plot/plot_data_17_8.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle