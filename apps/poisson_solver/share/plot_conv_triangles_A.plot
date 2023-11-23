set term postscript enhanced color
set output 'conv_triangles_A.eps'

set logscale x
set logscale y
set xrange [0.3:0.01]
set key bottom left
set grid xtics ytics mxtics

set title "Discrete energy error, triangles"

plot '../../../build/apps/poisson_solver/conv_triangles_0_std.txt'         using 1:3 w lp pt 4 title "k=0 standard", \
     '../../../build/apps/poisson_solver/conv_triangles_0_stabfree.txt'    using 1:3 w lp pt 5 title "k=0 stabfree", \
     '../../../build/apps/poisson_solver/conv_triangles_1_std.txt'         using 1:3 w lp pt 6 title "k=1 standard", \
     '../../../build/apps/poisson_solver/conv_triangles_1_stabfree.txt'    using 1:3 w lp pt 7 title "k=1 stabfree", \
     '../../../build/apps/poisson_solver/conv_triangles_2_std.txt'         using 1:3 w lp pt 8 title "k=2 standard", \
     '../../../build/apps/poisson_solver/conv_triangles_2_stabfree.txt'    using 1:3 w lp pt 9 title "k=2 stabfree", \
     '../../../build/apps/poisson_solver/conv_triangles_3_std.txt'         using 1:3 w lp pt 12 title "k=3 standard", \
     '../../../build/apps/poisson_solver/conv_triangles_3_stabfree.txt'    using 1:3 w lp pt 13 title "k=3 stabfree"
