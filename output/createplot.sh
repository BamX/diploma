#!/usr/local/bin/gnuplot -persist
set encoding utf8
set term png size 1200,700
set output 'plot.png'

set grid
set title 'График температуры слитка'
set xlabel 'Время, с'
set ylabel 'Температура, K'

set datafile separator ','
plot 'view.csv' using 1:2 title '(0.0, 0.115)' with lines, \
     'view.csv' using 1:3 title '(0.14, 0.115)' with lines, \
     'view.csv' using 1:4 title '(0.14, 0.0)' with lines, \
     'view.csv' using 1:5 title '(0.0, 0.0)' with lines, \
     'view.csv' using 1:6 title '(0.15, 0.125)' with lines
