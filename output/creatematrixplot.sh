#!/usr/local/bin/gnuplot -persist
set encoding utf8
set terminal gif animate delay 20 size 900,800
set output 'matrix.gif'
stats 'matrix.csv' nooutput

set xrange [0:100]
set yrange [0:100]
set cbrange [700:1800]

do for [i=1:int(STATS_blank)] {
    set title sprintf('Слиток time=%.2f с', (i-1) * 1816.0 / int(STATS_blank))
    plot 'matrix.csv' index (i-1) matrix with image
}
