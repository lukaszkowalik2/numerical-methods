set terminal svg
set output 'timing_comparison_large.svg'
set title 'Czas wykonania Sherman-Morrison dla dużych macierzy'
set xlabel 'Rozmiar macierzy'
set ylabel 'Czas wykonania [s]'
set grid
plot 'timing_data_large.txt' using 1:2 with linespoints title 'Sherman-Morrison'
