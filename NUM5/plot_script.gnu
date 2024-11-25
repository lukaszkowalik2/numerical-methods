set terminal svg
set output 'convergence_d3.svg'
set title 'Błąd w kolejnych iteracjach (d = 1.201)'
set xlabel 'Iteracja'
set ylabel 'log10(Błąd)'
plot 'errors_d3.dat' using 1:2 title 'Jacobi' with lines, \
     'errors_d3.dat' using 1:3 title 'Gauss-Seidel' with lines
