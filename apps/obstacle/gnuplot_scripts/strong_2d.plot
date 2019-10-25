set term postscript enhanced color font 'Helvetica,22'
set output 'strong_2d.eps'

set logscale x
set logscale y

set grid xtics mxtics ytics

set xrange [1:0.01] reverse

set style line 1 lt 1 lc rgb "#006600" lw 1
set style line 2 lt 3 dt 2 lc rgb "#990000" lw 1

plot 'triangles_k0.txt' using 1:2 w lp ls 1  pt 8   ps 1.5  title "Triangles, k=0", \
     'quads_k0.txt'     using 1:2 w lp ls 1  pt 4   ps 1.5  title "Quads, k=0", \
     'hexagons_k0.txt'  using 1:2 w lp ls 1  pt 14  ps 1.5  title "Hexagons, k=0", \
     'triangles_k1.txt' using 1:2 w lp ls 2  pt 8   ps 1.5  title "Triangles, k=1", \
     'quads_k1.txt'     using 1:2 w lp ls 2  pt 4   ps 1.5  title "Quads, k=1", \
     'hexagons_k1.txt'  using 1:2 w lp ls 2  pt 14  ps 1.5  title "Hexagons, k=1"
