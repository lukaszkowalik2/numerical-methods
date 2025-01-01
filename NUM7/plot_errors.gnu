set terminal svg enhanced size 1200,800
set grid
set key outside right spacing 2
set logscale y
set format y '10^{%L}'
set xrange [-1:1]
set yrange [1e-5:1e0]
set lmargin 12
set rmargin 6
set ylabel offset -2,0
set ytics offset 0,0

set output 'images/error_runge_n10.svg'
set title 'Błędy interpolacji funkcji f(x) = 1/(1+10x^2)\n = 10'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_runge_n10.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_runge_n20.svg'
set title 'Błędy interpolacji funkcji f(x) = 1/(1+10x^2)\n = 20'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_runge_n20.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_runge_n30.svg'
set title 'Błędy interpolacji funkcji f(x) = 1/(1+10x^2)\n = 30'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_runge_n30.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_sin3x_n10.svg'
set title 'Błędy interpolacji funkcji f(x) = sin(3x) + x\n = 10'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_sin3x_n10.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_sin3x_n20.svg'
set title 'Błędy interpolacji funkcji f(x) = sin(3x) + x\n = 20'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_sin3x_n20.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_sin3x_n30.svg'
set title 'Błędy interpolacji funkcji f(x) = sin(3x) + x\n = 30'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_sin3x_n30.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_sin_composite_n10.svg'
set title 'Błędy interpolacji funkcji f(x) = sin(x) + 0.3sin(50x)\n = 10'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_sin_composite_n10.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_sin_composite_n20.svg'
set title 'Błędy interpolacji funkcji f(x) = sin(x) + 0.3sin(50x)\n = 20'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_sin_composite_n20.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_sin_composite_n30.svg'
set title 'Błędy interpolacji funkcji f(x) = sin(x) + 0.3sin(50x)\n = 30'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_sin_composite_n30.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_rational_sin_n10.svg'
set title 'Błędy interpolacji funkcji f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2)\n = 10'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_rational_sin_n10.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_rational_sin_n20.svg'
set title 'Błędy interpolacji funkcji f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2)\n = 20'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_rational_sin_n20.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_rational_sin_n30.svg'
set title 'Błędy interpolacji funkcji f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2)\n = 30'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_rational_sin_n30.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_gauss_sin_n10.svg'
set title 'Błędy interpolacji funkcji f(x) = e^{-x^2}sin(10x)\n = 10'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_gauss_sin_n10.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_gauss_sin_n20.svg'
set title 'Błędy interpolacji funkcji f(x) = e^{-x^2}sin(10x)\n = 20'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_gauss_sin_n20.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

set output 'images/error_gauss_sin_n30.svg'
set title 'Błędy interpolacji funkcji f(x) = e^{-x^2}sin(10x)\n = 30'
set xlabel 'x'
set ylabel 'Błąd' offset -3
plot 'data/interpolation_gauss_sin_n30.dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'

unset logscale y
