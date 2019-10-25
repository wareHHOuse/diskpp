set term postscript enhanced color font 'Helvetica,22'
set output 'strong_3d.eps'

set logscale x
set logscale y

set grid xtics mxtics ytics

set xrange [1:0.01] reverse

set style line 1 lt 1 lc rgb "#006600" lw 1
set style line 2 lt 3 dt 2 lc rgb "#990000" lw 1

plot 'tetras_k0.txt' using 1:2 w lp ls 1  pt 8   ps 1.5  title "Tetrahedrons, k=0", \
     'hexas_k0.txt'  using 1:2 w lp ls 1  pt 4   ps 1.5  title "Cubes, k=0", \
     'prisms_k0.txt' using 1:2 w lp ls 1  pt 14  ps 1.5  title "Prisms, k=0", \
     'tetras_k1.txt' using 1:2 w lp ls 2  pt 8   ps 1.5  title "Tetrahedrons, k=1", \
     'hexas_k1.txt'  using 1:2 w lp ls 2  pt 4   ps 1.5  title "Cubes, k=1", \
     'prisms_k1.txt' using 1:2 w lp ls 2  pt 14  ps 1.5  title "Prisms, k=1"
