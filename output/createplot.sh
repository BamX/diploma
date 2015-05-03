#!/usr/local/bin/gnuplot -persist
set encoding utf8
set term png size 800,400
set output 'plot.png'

set grid
set title 'График охлаждения слитка'
set xlabel 'Время, с'
set ylabel 'Температура, с'

set datafile separator ','
plot 'view.csv'
