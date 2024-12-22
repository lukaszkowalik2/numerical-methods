set terminal svg enhanced size 800,600
set output 'qr_comparison.svg'
set title 'Convergence Comparison: Basic QR vs Wilkinson Shift'
set xlabel 'Iteration'
set ylabel 'log10(Error)'
plot 'qr_comparison.dat' using 1:2 title 'Basic QR' with lines lw 2,\
     'qr_comparison.dat' using 1:3 title 'Wilkinson Shift' with lines lw 2
