set terminal svg enhanced size 800,600
set output 'power_method_convergence.svg'
set title 'Power Method Convergence'
set xlabel 'Iteration'
set ylabel 'Error (log scale)'
set grid
plot 'power_method_convergence.dat' with lines lw 2 lt rgb 'blue' notitle