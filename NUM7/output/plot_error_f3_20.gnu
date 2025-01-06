set terminal svg enhanced size 800,600
set output 'errors_f3_20.svg'
set title 'Interpolation Errors for f3 (N=20)'
set xlabel 'x'
set ylabel 'Error'
set grid
set logscale y
plot 'errors_f3_20.dat' using 1:2 w l title 'Lagrange Error', \
     'errors_f3_20.dat' using 1:3 w l title 'Spline Error'
