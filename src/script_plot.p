set   autoscale
plot \
  'plot.dat' using 1:2       with lines lc rgb "black" lw 2 notitle,\
  'plot.dat' using 1:2:(0.6) with circles fill solid lc rgb "red" notitle,\
  'plot.dat' using 1:2:3     with labels tc rgb "white" offset (0,0) font 'Arial Bold' notitle
