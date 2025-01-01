set terminal svg enhanced size 1200,800
set grid
set key outside right spacing 2

set output 'images/interpolation_runge_n10.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = 1/(1+10x^2), n = 10 węzłów'
plot 'data/interpolation_runge_n10.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_runge_n20.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = 1/(1+10x^2), n = 20 węzłów'
plot 'data/interpolation_runge_n20.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_runge_n30.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = 1/(1+10x^2), n = 30 węzłów'
plot 'data/interpolation_runge_n30.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_sin3x_n10.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = sin(3x) + x, n = 10 węzłów'
plot 'data/interpolation_sin3x_n10.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_sin3x_n20.svg'
set xrange [-2:2]
set yrange [-2:2]
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = sin(3x) + x, n = 20 węzłów'
plot 'data/interpolation_sin3x_n20.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_sin3x_n30.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = sin(3x) + x, n = 30 węzłów'
plot 'data/interpolation_sin3x_n30.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_sin_composite_n10.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = sin(x) + 0.3sin(50x), n = 10 węzłów'
plot 'data/interpolation_sin_composite_n10.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_sin_composite_n20.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = sin(x) + 0.3sin(50x), n = 20 węzłów'
plot 'data/interpolation_sin_composite_n20.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_sin_composite_n30.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = sin(x) + 0.3sin(50x), n = 30 węzłów'
plot 'data/interpolation_sin_composite_n30.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_rational_sin_n10.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2), n = 10 węzłów'
plot 'data/interpolation_rational_sin_n10.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_rational_sin_n20.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2), n = 20 węzłów'
plot 'data/interpolation_rational_sin_n20.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_rational_sin_n30.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2), n = 30 węzłów'
plot 'data/interpolation_rational_sin_n30.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_gauss_sin_n10.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = e^{-x^2}sin(10x), n = 10 węzłów'
plot 'data/interpolation_gauss_sin_n10.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_gauss_sin_n20.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = e^{-x^2}sin(10x), n = 20 węzłów'
plot 'data/interpolation_gauss_sin_n20.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

set output 'images/interpolation_gauss_sin_n30.svg'
set xrange [-1:1]
set autoscale y
set xlabel 'x'
set ylabel 'y'
set title 'f(x) = e^{-x^2}sin(10x), n = 30 węzłów'
plot 'data/interpolation_gauss_sin_n30.dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \
     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \
     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'

