set term cairolatex color

set logscale x
set logscale y
set xrange [1:0.01]
set key bottom left
set format y "%1.1e"
set grid

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k0_eo_L2.tex"
set title "Triangles equal order k = 0 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_eo_0.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k0_moh_L2.tex"
set title "Triangles mixed order high k = 0 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_moh_0.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k1_eo_L2.tex"
set title "Triangles equal order k = 1 (L2)"
set yrange[1e-6:1e3]

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_eo_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k1_moh_L2.tex"
set title "Triangles mixed order high k = 1 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_moh_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k1_mol_L2.tex"
set title "Triangles mixed order low k = 1 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_mol_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k2_eo_L2.tex"
set title "Triangles equal order k = 2 (L2)"
set yrange [1e-9:1e1]

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_eo_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k2_moh_L2.tex"
set title "Triangles mixed order high k = 2 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_moh_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_tri_k2_mol_L2.tex"
set title "Triangles mixed order low k = 2 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_triangles_mol_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_triangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_triangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_triangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_triangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_triangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
########################################
########################################
set output "aniso_nostab/gnuplot/aniso_quad_k0_eo_L2.tex"
set title "Quads equal order k = 0 (L2)"
set yrange [1e-4:1]

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_eo_0.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_eo_0.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k0_moh_L2.tex"
set title "Quads mixed order high k = 0 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_moh_0.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_moh_0.txt'  using 1:2 w lp title "kyy = 0.00001" 

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k1_eo_L2.tex"
set title "Quads equal order k = 1 (L2)"
set yrange [1e-7:1e-1]

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_eo_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_eo_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k1_moh_L2.tex"
set title "Quads mixed order high k = 1 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_moh_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_moh_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k1_mol_L2.tex"
set title "Quads mixed order low k = 1 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_mol_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_mol_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k2_eo_L2.tex"
set title "Quads equal order k = 2 (L2)"
set yrange [1e-9:1e-2]

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_eo_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_eo_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k2_moh_L2.tex"
set title "Quads mixed order high k = 2 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_moh_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_moh_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_quad_k2_mol_L2.tex"
set title "Quads mixed order low k = 2 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_quadrangles_mol_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_quadrangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_quadrangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_quadrangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_quadrangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_quadrangles_mol_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
########################################
########################################
set output "aniso_nostab/gnuplot/aniso_hex_k0_eo_L2.tex"
set title "Hexagons equal order k = 0 (L2)"
set yrange [1e-5:1e2]

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_eo_0.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_eo_0.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_eo_0.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_eo_0.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_eo_0.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_eo_0.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k0_moh_L2.tex"
set title "Hexagons mixed order high k = 0 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_moh_0.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_moh_0.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_moh_0.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_moh_0.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_moh_0.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_moh_0.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k1_eo_L2.tex"
set title "Hexagons equal order k = 1 (L2)"
set yrange [1e-7:1e1]

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_eo_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_eo_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_eo_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_eo_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_eo_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_eo_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k1_moh_L2.tex"
set title "Hexagons mixed order high k = 1 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_moh_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_moh_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_moh_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_moh_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_moh_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_moh_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k1_mol_L2.tex"
set title "Hexagons mixed order low k = 1 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_mol_1.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_mol_1.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_mol_1.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_mol_1.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_mol_1.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_mol_1.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k2_eo_L2.tex"
set title "Hexagons equal order k = 2 (L2)"
set yrange [1e-9:1e-2]

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_eo_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_eo_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_eo_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_eo_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_eo_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_eo_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k2_moh_L2.tex"
set title "Hexagons mixed order high k = 2 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_moh_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_moh_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_moh_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_moh_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_moh_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_moh_2.txt'  using 1:2 w lp title "kyy = 0.00001"

########################################
set output "aniso_nostab/gnuplot/aniso_hex_k2_mol_L2.tex"
set title "Hexagons mixed order low k = 2 (L2)"

plot 'aniso_nostab/data/aniso_1.0_conv_hexagons_mol_2.txt'  using 1:2 w lp title "kyy = 1.0",\
     'aniso_nostab/data/aniso_0.1_conv_hexagons_mol_2.txt'  using 1:2 w lp title "kyy = 0.1",\
     'aniso_nostab/data/aniso_0.01_conv_hexagons_mol_2.txt'  using 1:2 w lp title "kyy = 0.01",\
     'aniso_nostab/data/aniso_0.001_conv_hexagons_mol_2.txt'  using 1:2 w lp title "kyy = 0.001",\
     'aniso_nostab/data/aniso_0.0001_conv_hexagons_mol_2.txt'  using 1:2 w lp title "kyy = 0.0001",\
     'aniso_nostab/data/aniso_1e-05_conv_hexagons_mol_2.txt'  using 1:2 w lp title "kyy = 0.00001"

