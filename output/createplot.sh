#!/usr/local/bin/gnuplot -persist
set encoding utf8
set term png size 800,400
set output 'plot.png'

set grid
set title 'График охлаждения слитка'
set xlabel 'Время, с'
set ylabel 'Температура, с'

set datafile separator ','
plot 'view.csv' using 1:2 title 'центр, низ' with lines, \
     'view.csv' using 1:3 title 'центр, центр' with lines, \
     'view.csv' using 1:4 title 'лево, центр' with lines, \
     'view.csv' using 1:5 title 'ноль' with lines, \
     'view.csv' using 1:6 title 'граница' with lines
