set terminal svg enhanced size 800,600
set output 'interpolation_y_10.svg'
set title 'Interpolation for y (N=10)'
set xlabel 'x'
set ylabel 'y'
set grid
plot 'exact_y_10.dat' w l title 'Exact', \
     'lagrange_y_10.dat' w l title 'Lagrange', \
     'spline_y_10.dat' w l title 'Spline', \
     'points_y_10.dat' w p pt 7 title 'Points'
