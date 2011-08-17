#!/bin/sh
gnuplot << EOF
set terminal postscript eps color enhanced
set output "$1.eps"
set xlabel "Energy [MeV]"
set ylabel "Cross Section [b]"
set title "(n,2n) reaction"
set xrange [ 0 : 20 ]
set yrange [ 0 : 2 ]
set mxtics 5
set mytics 5
set xtics 5
set ytics 0.5
plot "$1.dat" using 1:2 notitle w l
EOF


set yrange [0:*]
set tics scale 0.0
set ytics







set size 1, 1
set term png size 800, 600
set output "12-10F.20110119B_hiseq2000.PE.filtered.summary.pct_gc_histo.png"

set title "Percent GC of sequencing reads\n(histogram)"
set boxwidth 0.9 relative
set style fill solid 1.0 border -1
set style line 1 lt 2 lc rgb "blue" lw 2 pt 2 ps 0.5

set key outside below
set key box lw 1 lc rgb "black"
set key spacing 1.5

set xtics nomirror 
set xlabel "Percent GC"

set ytics nomirror 
set y2tics nomirror 

set yrange [0:*]
set ytics border
set mytics
set ylabel "Counts"

set y2range [0:*]
set y2tics border
set y2tics
set my2tics
set y2label "Percent"

plot "12-10F.20110119B_hiseq2000.PE.filtered.summary.pct_gc_histo.txt" using 1:2 every ::2 title "Read counts" with boxes, "12-10F.20110119B_hiseq2000.PE.filtered.summary.pct_gc_histo.txt" using 1:3 every ::2 axes x1y2 title "Percentage of reads" with linespoints ls 1




set size 1, 1
set term png size 800, 600
set output "12-10F.merged.recal.realn.reSrt.coverage.coverage_histogram_per_chromosome.man.png"

set title "Percent GC of sequencing reads\n(histogram)"
set boxwidth 0.9 relative
set style fill solid 1.0 border -1
set style line 1 lt 2 lc rgb "blue" lw 2 pt 2 ps 0.5

set key outside below
set key box lw 1 lc rgb "black"
set key spacing 1.5

set xtics nomirror 
set xlabel "Percent GC"

set ytics nomirror 
set y2tics nomirror 

set yrange [0:*]
set ytics border
set mytics
set ylabel "Counts"

set y2range [0:*]
set y2tics border
set y2tics
set my2tics
set y2label "Percent"

set key autotitle columnhead

plot "12-10F.merged.recal.realn.reSrt.coverage.coverage_histogram_per_chromosome.txt" with linespoints














