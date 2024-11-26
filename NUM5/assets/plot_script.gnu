set terminal svg size 800,600
set output 'assets/convergence_d5.000000_start100.000000.svg'
set title 'Błąd w kolejnych iteracjach (d = 5)'
set xlabel 'Iteracja'
set ylabel 'log10(Błąd)'
set grid
plot 'assets/errors_d5.000000_start100.000000.dat' using 1:2 title 'Jacobi' with lines lw 2, \
     'assets/errors_d5.000000_start100.000000.dat' using 1:3 title 'Gauss-Seidel' with lines lw 2
