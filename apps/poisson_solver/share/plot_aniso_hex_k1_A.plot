set logscale x
set logscale y
set xrange [1:0.01]
set key bottom left
set title "Triangles k=1 equal order"

plot '../../../build/apps/poisson_solver/aniso_1.0_conv_hexagons_1.txt' using 1:3 w lp title "k_{yy} = 1", \
     '../../../build/apps/poisson_solver/aniso_0.1_conv_hexagons_1.txt' using 1:3 w lp title "k_{yy} = 0.1", \
     '../../../build/apps/poisson_solver/aniso_0.01_conv_hexagons_1.txt' using 1:3 w lp title "k_{yy} = 0.01", \
     '../../../build/apps/poisson_solver/aniso_0.001_conv_hexagons_1.txt' using 1:3 w lp title "k_{yy} = 0.001"
