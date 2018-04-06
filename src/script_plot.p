set   autoscale
plot \
  'plot.dat' using 1:2       with lines lc rgb "green" lw 2 notitle,\
  'plot.dat' using 1:2:(0.6) with circles fill solid lc rgb "black" notitle,\
  'plot.dat' using 1:2:3     with labels tc rgb "black" offset (0,0) font 'Arial Bold' notitle
