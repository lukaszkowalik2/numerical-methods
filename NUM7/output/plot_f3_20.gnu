set terminal svg enhanced size 800,600
set output 'interpolation_f3_20.svg'
set title 'Interpolation for f3 (N=20)'
set xlabel 'x'
set ylabel 'y'
set grid
plot 'exact_f3_20.dat' w l title 'Exact', \
     'lagrange_f3_20.dat' w l title 'Lagrange', \
     'spline_f3_20.dat' w l title 'Spline', \
     'points_f3_20.dat' w p pt 7 title 'Points'
