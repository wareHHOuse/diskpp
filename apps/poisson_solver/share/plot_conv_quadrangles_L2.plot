set term postscript enhanced color
set output 'conv_quadrangles_L2.eps'

set logscale x
set logscale y
set xrange [0.3:0.01]
set key bottom left

plot '../../../build/apps/poisson_solver/conv_quadrangles_0_std.txt' using 1:2 w lp title "k=0 standard", \
     '../../../build/apps/poisson_solver/conv_quadrangles_0_stabfree.txt' using 1:2 w lp title "k=0 stabfree", \
     '../../../build/apps/poisson_solver/conv_quadrangles_1_std.txt' using 1:2 w lp title "k=1 standard", \
     '../../../build/apps/poisson_solver/conv_quadrangles_1_stabfree.txt' using 1:2 w lp title "k=1 stabfree", \
     '../../../build/apps/poisson_solver/conv_quadrangles_2_std.txt' using 1:2 w lp title "k=2 standard", \
     '../../../build/apps/poisson_solver/conv_quadrangles_2_stabfree.txt' using 1:2 w lp title "k=2 stabfree", \
     '../../../build/apps/poisson_solver/conv_quadrangles_3_std.txt' using 1:2 w lp title "k=3 standard", \
     '../../../build/apps/poisson_solver/conv_quadrangles_3_stabfree.txt' using 1:2 w lp title "k=3 stabfree"
