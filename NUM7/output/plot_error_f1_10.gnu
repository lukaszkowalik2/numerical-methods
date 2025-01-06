set terminal svg enhanced size 800,600
set output 'errors_f1_10.svg'
set title 'Interpolation Errors for f1 (N=10)'
set xlabel 'x'
set ylabel 'Error'
set grid
set logscale y
plot 'errors_f1_10.dat' using 1:2 w l title 'Lagrange Error', \
     'errors_f1_10.dat' using 1:3 w l title 'Spline Error'
