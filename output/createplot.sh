#!/usr/local/bin/gnuplot -persist
set encoding utf8
set term png size 800,400
set output 'plot.png'

set grid
set title 'График охлаждения слитка'
set xlabel 'Время, с'
set ylabel 'Температура, с'

set datafile separator ','
plot 'view.csv' using 1:2 title '5 мм от середины' with lines, \
     'view.csv' using 1:3 title 'между центром и серединой' with lines, \
     'view.csv' using 1:4 title '5 мм от поверхности' with lines
