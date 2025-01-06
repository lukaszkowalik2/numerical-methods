set terminal svg enhanced size 800,600
set output 'interpolation_y_30.svg'
set title 'Interpolation for y (N=30)'
set xlabel 'x'
set ylabel 'y'
set grid
set yrange [-2:2]
plot 'exact_y_30.dat' w l title 'Exact', \
     'lagrange_y_30.dat' w l title 'Lagrange', \
     'spline_y_30.dat' w l title 'Spline', \
     'points_y_30.dat' w p pt 7 title 'Points'
