set terminal svg
set output 'timing_comparison.svg'
set title 'Porównanie czasów wykonania'
set xlabel 'Rozmiar macierzy'
set ylabel 'Czas wykonania [s]'
set grid
plot 'timing_data.txt' using 1:2 with linespoints title 'Sherman-Morrison', 'timing_data.txt' using 1:3 with linespoints title 'Eigen', 'timing_data.txt' using 1:4 with linespoints title 'Gauss', 'timing_data.txt' using 1:5 with linespoints title 'QR'
