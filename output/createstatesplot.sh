#!/usr/local/bin/gnuplot -persist
set encoding utf8
set terminal gif animate delay 5 size 800,400
set output 'states.gif'
stats 'states.out' nooutput

set xrange [0:100]
set yrange [0:106]
set cbrange [0:2]

do for [i=1:int(STATS_blank)] {
    set title sprintf('%d', (i-1))
    plot 'states.out' index (i-1) matrix with image
}
