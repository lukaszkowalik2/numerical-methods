set terminal svg enhanced size 1200,800
set output 'diagonal_evolution.svg'
set title 'Zbieznosc do dokladnych wartosci wlasnych'
set xlabel 'Iteracja'
set ylabel 'log10(|λi - λi^{exact}|)'
set xrange [-1:65]
set grid
plot 'diagonal_evolution.dat' using 1:2 title '|λ1 - λ1^{exact}|' with lines lw 2,\
     'diagonal_evolution.dat' using 1:3 title '|λ2 - λ2^{exact}|' with lines lw 2,\
     'diagonal_evolution.dat' using 1:4 title '|λ3 - λ3^{exact}|' with lines lw 2,\
     'diagonal_evolution.dat' using 1:5 title '|λ4 - λ4^{exact}|' with lines lw 2
