set terminal pngcairo enhanced size 800,600
set output 'plot.png'

set title 'Aproksymacja'
set xlabel 'x'
set ylabel 'y'

set grid

# Plot points, original function and approximation
plot 'original.dat' using 1:2 with lines lw 2 title 'Funkcja F(x)', \
     'points.dat' using 1:2 with points pt 7 ps 0.5 title 'Punkty pomiarowe', \
     'approximation.dat' using 1:2 with lines lw 2 dt 2 title 'Aproksymacja'
